/**********************************************************************************

 Copyright (c) 2020 Jonas Sauer

 MIT License

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
 modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 is furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
 IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

**********************************************************************************/

#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Helpers/Vector/Vector.h"

#include "InitialTransfers.h"
#include "Profiler.h"

#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/Container/Set.h"
#include "../../DataStructures/Container/Map.h"
#include "../../DataStructures/Container/ExternalKHeap.h"

namespace RAPTOR {

template<typename INITIAL_TRANSFERS, typename PROFILER>
class OneToAllMCR {

public:
    using InitialTransferType = INITIAL_TRANSFERS;
    using InitialTransferGraph = typename InitialTransferType::Graph;
    using Profiler = PROFILER;
    using Type = OneToAllMCR<InitialTransferType, Profiler>;

public:
    struct EarliestArrivalLabel {
        EarliestArrivalLabel() : arrivalTime(never), parentDepartureTime(never), parent(noVertex), usesRoute(false), routeId(noRouteId), timestamp(0) {}
        int arrivalTime;
        int parentDepartureTime;
        Vertex parent;
        bool usesRoute;
        RouteId routeId;
        int timestamp;
    };
    using Round = std::vector<EarliestArrivalLabel>;

    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel() : arrivalTime(never), parent(noVertex), timestamp(0) {}
        int arrivalTime;
        Vertex parent;
        int timestamp;
        inline bool hasSmallerKey(const DijkstraLabel* const other) const noexcept {
            return arrivalTime < other->arrivalTime;
        }
    };

public:
    template<typename ATTRIBUTE>
    OneToAllMCR(const Data& data, const InitialTransferGraph& forwardGraph, const InitialTransferGraph& backwardGraph, const ATTRIBUTE weight, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialTransfers(forwardGraph, backwardGraph, data.numberOfStops(), weight),
        roundIndex(-1),
        stopsUpdatedByRoute(data.numberOfStops()),
        stopsUpdatedByTransfer(data.numberOfStops()),
        routesServingUpdatedStops(data.numberOfRoutes()),
        sourceVertex(noVertex),
        targetVertex(noVertex),
        targetStop(noStop),
        sourceDepartureTime(intMax),
        timestamp(0),
        labelByNumberOfTrips(1, std::vector<DijkstraLabel>(data.transferGraph.numVertices())),
        profiler(profilerTemplate) {
        AssertMsg(data.hasImplicitBufferTimes(), "Either min transfer times have to be used OR departure buffer times have to be implicit!");
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_VERTICES, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    template<typename T = CHGraph, typename = std::enable_if_t<Meta::Equals<T, CHGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    OneToAllMCR(const Data& data, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) :
        OneToAllMCR(data, chData.forward, chData.backward, Weight, profilerTemplate) {
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    OneToAllMCR(const Data& data, const TransferGraph& forwardGraph, const TransferGraph& backwardGraph, const Profiler& profilerTemplate = Profiler()) :
        OneToAllMCR(data, forwardGraph, backwardGraph, TravelTime, profilerTemplate) {
    }

    template<bool CLEAR = true>
    inline void run(const Vertex source, const int departureTime, const Vertex target = noVertex, const size_t maxRounds = 50) noexcept {
        runInitialize<CLEAR>(source, departureTime, target);
        runInitialTransfers();
        evaluateInitialTransfers();
        runRounds(maxRounds);
    }

    template<bool CLEAR = true>
    inline void runInitialize(const Vertex source, const int departureTime, const Vertex target = noVertex) noexcept {
        profiler.start();
        profiler.startExtraRound(EXTRA_ROUND_CLEAR);
        clear<CLEAR>();
        profiler.doneRound();

        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        initialize(source, departureTime, target);
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.doneRound();
    }

    inline void runInitialTransfers() noexcept {
        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        relaxInitialTransfers();
        profiler.donePhase(PHASE_TRANSFERS);
        profiler.doneRound();
    }

    inline void runAddSource(const StopId source, const int departureTime) noexcept {
        EarliestArrivalLabel& label = currentRoundLabel(source);
        label.arrivalTime = departureTime;
        label.parentDepartureTime = sourceDepartureTime;
        label.usesRoute = false;
        stopsUpdatedByTransfer.insert(source);
        dijkstraLabel(roundIndex, source).arrivalTime = departureTime;
    }

    inline void runRounds(const size_t maxRounds = 50) noexcept {
        profiler.startRound();
        profiler.startPhase();
        startNewRound();
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.startPhase();
        collectRoutesServingUpdatedStops();
        profiler.donePhase(PHASE_COLLECT);
        profiler.startPhase();
        scanRoutes();
        profiler.donePhase(PHASE_SCAN);
        for (size_t i = 0; (i < maxRounds) && (!stopsUpdatedByRoute.empty()); i++) {
            profiler.startPhase();
            relaxIntermediateTransfers();
            profiler.donePhase(PHASE_TRANSFERS);
            profiler.doneRound();
            profiler.startRound();
            profiler.startPhase();
            startNewRound();
            profiler.donePhase(PHASE_INITIALIZATION);
            profiler.startPhase();
            collectRoutesServingUpdatedStops();
            profiler.donePhase(PHASE_COLLECT);
            profiler.startPhase();
            scanRoutes();
            profiler.donePhase(PHASE_SCAN);
        }
        profiler.doneRound();
        profiler.done();
    }

    template<bool RESET = true>
    inline void clear() noexcept {
        roundIndex = -1;
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        sourceVertex = noVertex;
        targetVertex = noVertex;
        targetStop = StopId(data.numberOfStops());
        queue.clear();
        if constexpr (RESET) {
            timestamp = 1;
            std::vector<Round>().swap(rounds);
            std::vector<std::vector<DijkstraLabel>>(1).swap(labelByNumberOfTrips);
            std::vector<DijkstraLabel>(data.transferGraph.numVertices()).swap(labelByNumberOfTrips[0]);
        } else {
            timestamp++;
        }
    }

    inline void reset() noexcept {
        clear<true>();
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void initialize(const Vertex source, const int departureTime, const Vertex target) noexcept {
        sourceVertex = source;
        targetVertex = target;
        sourceDepartureTime = departureTime;
        startNewRound();
        if (data.isStop(source)) {
            stopsUpdatedByRoute.insert(StopId(source));
            EarliestArrivalLabel& label = currentRoundLabel(StopId(source));
            label.arrivalTime = departureTime;
            label.parent = source;
            label.parentDepartureTime = departureTime;
            label.usesRoute = false;
        }
    }

    inline void evaluateInitialTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        for (const Vertex stop : initialTransfers.getForwardPOIs()) {
            AssertMsg(data.isStop(stop), "Reached POI " << stop << " is not a stop!");
            AssertMsg(initialTransfers.getForwardDistance(stop) != INFTY, "Vertex " << stop << " was not reached!");
            //The initial transfers are evaluated automatically when the label is updated
            roundLabel(0, StopId(stop));
            stopsUpdatedByTransfer.insert(StopId(stop));
        }
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        for (const StopId stop : stopsUpdatedByTransfer) {
            for (const RouteSegment& route : data.routesContainingStop(stop)) {
                AssertMsg(data.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                AssertMsg(data.stopIds[data.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (route.stopIndex + 1 == data.numberOfStopsInRoute(route.routeId)) continue;
                if (data.lastTripOfRoute(route.routeId)[route.stopIndex].departureTime < previousRoundLabel(stop).arrivalTime) continue;
                if (routesServingUpdatedStops.contains(route.routeId)) {
                    routesServingUpdatedStops[route.routeId] = std::min(routesServingUpdatedStops[route.routeId], route.stopIndex);
                } else {
                    routesServingUpdatedStops.insert(route.routeId, route.stopIndex);
                }
            }
        }
    }

    inline void scanRoutes() noexcept {
        stopsUpdatedByRoute.clear();
        for (const RouteId route : routesServingUpdatedStops.getKeys()) {
            profiler.countMetric(METRIC_ROUTES);
            StopIndex stopIndex = routesServingUpdatedStops[route];
            const size_t tripSize = data.numberOfStopsInRoute(route);
            AssertMsg(stopIndex < tripSize - 1, "Cannot scan a route starting at/after the last stop (Route: " << route << ", StopIndex: " << stopIndex << ", TripSize: " << tripSize << ", RoundIndex: " << roundIndex << ")!");

            const StopId* stops = data.stopArrayOfRoute(route);
            const StopEvent* firstTrip = data.firstTripOfRoute(route);
            const StopEvent* trip = data.lastTripOfRoute(route);
            StopId stop = stops[stopIndex];
            AssertMsg(trip[stopIndex].departureTime >= previousRoundLabel(stop).arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << stop << ", StopIndex: " << stopIndex << ", Time: " << previousRoundLabel(stop).arrivalTime << ", LastDeparture: " << trip[stopIndex].departureTime << ", RoundIndex: " << roundIndex << ")!");

            StopIndex parentIndex = stopIndex;
            while (stopIndex < tripSize - 1) {
                while ((trip > firstTrip) && ((trip - tripSize + stopIndex)->departureTime >= previousRoundLabel(stop).arrivalTime)) {
                    trip -= tripSize;
                    parentIndex = stopIndex;
                }
                stopIndex++;
                stop = stops[stopIndex];
                profiler.countMetric(METRIC_ROUTE_SEGMENTS);
                if (arrivalByRoute(stop, trip[stopIndex].arrivalTime)) {
                    EarliestArrivalLabel& label = currentRoundLabel(stop);
                    label.parent = stops[parentIndex];
                    label.parentDepartureTime = trip[parentIndex].departureTime;
                    label.usesRoute = true;
                    label.routeId = route;
                }
            }
        }
    }

    inline void relaxInitialTransfers() noexcept {
        initialTransfers.template run<FORWARD, BACKWARD>(sourceVertex);
    }

    inline void relaxIntermediateTransfers() noexcept {
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        AssertMsg(queue.empty(), "Queue still has " << queue.size() << " elements!");
        for (const StopId stop : stopsUpdatedByRoute) {
            const int arrivalTime = currentRoundLabel(stop).arrivalTime;
            arrivalByEdge<true>(stop, arrivalTime, stop);
        }
        dijkstra();
    }

    inline void dijkstra() noexcept {
        while (!queue.empty()) {
            DijkstraLabel* uLabel = queue.extractFront();
            const Vertex u = Vertex(uLabel - &(labelByNumberOfTrips[roundIndex][0]));
            for (Edge edge : data.transferGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = data.transferGraph.get(ToVertex, edge);
                if (v == targetVertex || v == uLabel->parent) continue;
                arrivalByEdge<true>(v, uLabel->arrivalTime + data.transferGraph.get(TravelTime, edge), uLabel->parent);
            }
            if (data.isStop(u)) {
                arrivalByTransfer(StopId(u), uLabel->arrivalTime, uLabel->parent, getParentDepartureTime(uLabel->parent, roundIndex));
            }
            profiler.countMetric(METRIC_VERTICES);
        }
        if (!queue.empty()) {
            queue.clear();
        }
    }

    inline EarliestArrivalLabel& roundLabel(const size_t round, const StopId stop) noexcept {
        EarliestArrivalLabel& result = rounds[round][stop];
        if (result.timestamp != timestamp) {
            if (round >= 1) {
                result.arrivalTime = std::min(result.arrivalTime, roundLabel(round - 1, stop).arrivalTime);
            } else if (round == 0) {
                int distance = (stop == targetStop) ? initialTransfers.getDistance() : initialTransfers.getForwardDistance(stop);
                if (distance != INFTY) {
                    profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
                    result.arrivalTime = sourceDepartureTime + distance;
                    result.parent = sourceVertex;
                    result.parentDepartureTime = sourceDepartureTime;
                    result.usesRoute = false;
                }
            }
            result.timestamp = timestamp;
        }
        return result;
    }

    inline DijkstraLabel& dijkstraLabel(const size_t numTrips, const Vertex vertex) noexcept {
        DijkstraLabel& result = labelByNumberOfTrips[numTrips][vertex];
        if (result.timestamp != timestamp) {
            if (numTrips >= 1) {
                result.arrivalTime = std::min(result.arrivalTime, dijkstraLabel(numTrips - 1, vertex).arrivalTime);
            } else if (numTrips == 0) {
                int distance = (vertex == targetVertex) ? initialTransfers.getDistance() : initialTransfers.getForwardDistance(vertex);
                if (distance != INFTY) {
                    result.arrivalTime = sourceDepartureTime + distance;
                    result.parent = sourceVertex;
                }
            }
            result.timestamp = timestamp;
        }
        return result;
    }

    inline EarliestArrivalLabel& currentRoundLabel(const StopId stop) noexcept {
        AssertMsg(roundIndex < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        return roundLabel(roundIndex, stop);
    }

    inline EarliestArrivalLabel& previousRoundLabel(const StopId stop) noexcept {
        AssertMsg(roundIndex - 1 < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        AssertMsg(roundIndex > 0, "Cannot return previous round, because no round exists!");
        return roundLabel(roundIndex - 1, stop);
    }

    inline void startNewRound() noexcept {
        AssertMsg(roundIndex + 1 <= rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        roundIndex++;
        if (roundIndex == rounds.size()) {
            if (rounds.empty()) {
                rounds.emplace_back(data.numberOfStops());
            } else {
                rounds.emplace_back(rounds.back());
            }
        }
        if (roundIndex == labelByNumberOfTrips.size()) {
            labelByNumberOfTrips.emplace_back(labelByNumberOfTrips.back());
        }
    }

    inline bool arrivalByRoute(const StopId stop, const int arrivalTime) noexcept {
        AssertMsg(roundIndex > 0, "arrivalByRoute cannot be used in the first round!");
        AssertMsg(data.isStop(stop), "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], stop: " << stop << ")!");
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if ((label.arrivalTime <= arrivalTime) || (previousRoundLabel(stop).arrivalTime <= arrivalTime)) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        label.arrivalTime = arrivalTime;
        stopsUpdatedByRoute.insert(stop);
        return true;
    }

    template<bool ADD_TO_QUEUE>
    inline bool arrivalByEdge(const Vertex vertex, const int arrivalTime, const Vertex parent) noexcept {
        AssertMsg(data.isStop(parent), "Parent vertex (" << parent << ") is not a stop");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");
        DijkstraLabel& label = dijkstraLabel(roundIndex, vertex);
        if (label.arrivalTime <= arrivalTime) return false;

        label.arrivalTime = arrivalTime;
        label.parent = parent;
        if constexpr (ADD_TO_QUEUE) queue.update(&label);
        return true;
    }

    inline void arrivalByTransfer(const StopId stop, const int arrivalTime, const Vertex parent, const int parentDepartureTime) noexcept {
        AssertMsg(data.isStop(stop) || stop == targetStop, "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if (data.isStop(stop)) stopsUpdatedByTransfer.insert(stop);
        if (label.arrivalTime <= arrivalTime) return;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        label.arrivalTime = arrivalTime;
        label.parent = parent;
        label.parentDepartureTime = parentDepartureTime;
        label.usesRoute = false;
    }

    inline int getParentDepartureTime(const Vertex vertex, const size_t round) noexcept {
        AssertMsg(data.isStop(vertex) || vertex == sourceVertex, "Vertex " << vertex << " should be a stop or the source vertex!");
        if (vertex == sourceVertex) {
            return sourceDepartureTime;
        } else {
            return roundLabel(round, StopId(vertex)).arrivalTime;
        }
    }

private:
    const Data& data;
    TransferGraph minChangeTimeGraph;

    InitialTransferType initialTransfers;

    std::vector<Round> rounds;
    size_t roundIndex;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Vertex sourceVertex;
    Vertex targetVertex;
    StopId targetStop; //One-to-one only
    int sourceDepartureTime;
    int timestamp;

    std::vector<std::vector<DijkstraLabel>> labelByNumberOfTrips;
    ExternalKHeap<2, DijkstraLabel> queue;

    Profiler profiler;

};

}
