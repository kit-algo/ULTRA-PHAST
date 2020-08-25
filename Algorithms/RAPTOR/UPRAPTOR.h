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

#include "../CH/UPQuery/GroupedUPQuery.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/Container/Set.h"
#include "../../DataStructures/Container/Map.h"
#include "Profiler.h"

namespace RAPTOR {

template<size_t GROUPED_ROUNDS, bool USE_DFS_ORDER, typename PROFILER = NoProfiler>
class UPRAPTOR {

public:
    constexpr static size_t GroupedRounds = GROUPED_ROUNDS;
    constexpr static bool UseDFSOrder = USE_DFS_ORDER;
    using Profiler = PROFILER;
    constexpr static bool Debug = !Meta::Equals<Profiler, NoProfiler>();
    using Type = UPRAPTOR<GroupedRounds, UseDFSOrder, Profiler>;
    using InitialAndFinalTransfers = CH::GroupedUPQuery<true, Debug, GroupedRounds>;

private:
    struct EarliestArrivalLabel {
        EarliestArrivalLabel() : arrivalTime(never), parentDepartureTime(never), parent(noVertex), usesRoute(false), routeId(noRouteId) {}
        int arrivalTime;
        int parentDepartureTime;
        Vertex parent;
        bool usesRoute;
        union {
            RouteId routeId;
            Edge transferId;
        };
    };
    using Round = std::vector<EarliestArrivalLabel>;

    UPRAPTOR(const Data& data, const TransferGraph& transferGraph, const CH::CH& chData, const Order&& chOrder, const IndexedSet<false, Vertex>& targets, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialAndFinalTransfers(transferGraph, chData, std::move(chOrder), data.numberOfStops(), targets),
        sourceVertex(noVertex),
        sourceDepartureTime(never),
        targetVertices(targets),
        earliestArrival(data.numberOfStops()),
        stopsUpdatedByRoute(data.numberOfStops()),
        stopsUpdatedByTransfer(data.numberOfStops()),
        routesServingUpdatedStops(data.numberOfRoutes()),
        profiler(profilerTemplate) {
        AssertMsg(data.hasImplicitBufferTimes(), "Implicit departure buffer times must be used!");
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION, EXTRA_ROUND_FINAL_TRANSFERS});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS, PHASE_FINAL_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

public:
    inline static Order vertexOrder(const CH::CH& chData) noexcept {
        if constexpr (UseDFSOrder) {
            return Order(Vector::reverse(CH::getOrder(chData)));
        } else {
            return Order(CH::getLevelOrderTopDown(chData));
        }
    }

    inline static UPRAPTOR Reordered(Data& data, TransferGraph& transferGraph, CH::CH& chData, IndexedSet<false, Vertex>& targets, const Profiler& profilerTemplate = Profiler()) noexcept {
        Order chOrder = vertexOrder(chData);
        Order fullOrder = chOrder.splitAt(data.numberOfStops());
        Order stopOrder = fullOrder;
        stopOrder.resize(data.numberOfStops());

        data.applyStopOrder(stopOrder);
        transferGraph.applyVertexOrder(fullOrder);
        chData.applyVertexOrder(fullOrder);
        targets.applyPermutation(Permutation(Construct::Invert, fullOrder));

        size_t numStops = 0;
        size_t numVertices = data.numberOfStops();
        Order phastOrder;
        for (const size_t i : chOrder) {
            if (i < data.numberOfStops()) {
                phastOrder.emplace_back(numStops++);
            } else {
                phastOrder.emplace_back(numVertices++);
            }
        }

        return UPRAPTOR(data, transferGraph, chData, std::move(phastOrder), targets, profilerTemplate);
    }

    inline void run(const Vertex source, const int departureTime, const size_t maxRounds = 50) noexcept {
        profiler.start();
        profiler.startExtraRound(EXTRA_ROUND_CLEAR);
        clear();
        profiler.doneRound();

        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        initialize(source, departureTime);
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.startPhase();
        relaxInitialTransfers();
        profiler.donePhase(PHASE_TRANSFERS);
        profiler.doneRound();

        for (size_t i = 0; i < maxRounds; i++) {
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
            if (stopsUpdatedByRoute.empty()) {
                profiler.doneRound();
                break;
            }
            relaxIntermediateTransfers();
            profiler.doneRound();
        }
        profiler.startExtraRound(EXTRA_ROUND_FINAL_TRANSFERS);
        profiler.startPhase();
        initialAndFinalTransfers.finalize();
        profiler.donePhase(PHASE_FINAL_TRANSFERS);
        profiler.doneRound();
        profiler.done();
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

    inline long long getUpwardSweepGraphVertices() const noexcept {
        return initialAndFinalTransfers.getUpwardSweepGraphVertices();
    }

    inline long long getUpwardSweepGraphEdges() const noexcept {
        return initialAndFinalTransfers.getUpwardSweepGraphEdges();
    }

    inline long long getStopGraphVertices() const noexcept {
        return initialAndFinalTransfers.getStopGraphVertices();
    }

    inline long long getStopGraphEdges() const noexcept {
        return initialAndFinalTransfers.getStopGraphEdges();
    }

    inline long long getTargetGraphVertices() const noexcept {
        return initialAndFinalTransfers.getTargetGraphVertices();
    }

    inline long long getTargetGraphEdges() const noexcept {
        return initialAndFinalTransfers.getTargetGraphEdges();
    }

private:
    template<bool RESET_CAPACITIES = false>
    inline void clear() noexcept {
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        if constexpr (RESET_CAPACITIES) {
            std::vector<Round>().swap(rounds);
            std::vector<int>(earliestArrival.size(), never).swap(earliestArrival);
        } else {
            rounds.clear();
            Vector::fill(earliestArrival, never);
        }
    }

    inline void reset() noexcept {
        clear<true>();
    }

    inline void initialize(const Vertex source, const int departureTime) noexcept {
        sourceVertex = source;
        sourceDepartureTime = departureTime;
        initialAndFinalTransfers.initialize();
        startNewRound();
        if (data.isStop(source)) {
            arrivalByRoute(StopId(source), sourceDepartureTime);
            currentRound()[source].parent = source;
            currentRound()[source].parentDepartureTime = sourceDepartureTime;
            currentRound()[source].usesRoute = false;
            stopsUpdatedByTransfer.insert(StopId(source));
        }
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        for (const StopId stop : stopsUpdatedByTransfer) {
            AssertMsg(data.isStop(stop), "Stop " << stop << " is out of range!");
            const int arrivalTime = previousRound()[stop].arrivalTime;
            AssertMsg(arrivalTime < never, "Updated stop has arrival time = never!");
            for (const RouteSegment& route : data.routesContainingStop(stop)) {
                AssertMsg(data.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                AssertMsg(data.stopIds[data.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (route.stopIndex + 1 == data.numberOfStopsInRoute(route.routeId)) continue;
                if (data.lastTripOfRoute(route.routeId)[route.stopIndex].departureTime < arrivalTime) continue;
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
            AssertMsg(stopIndex < tripSize - 1, "Cannot scan a route starting at/after the last stop (Route: " << route << ", StopIndex: " << stopIndex << ", TripSize: " << tripSize << ")!");

            const StopId* stops = data.stopArrayOfRoute(route);
            const StopEvent* trip = data.lastTripOfRoute(route);
            StopId stop = stops[stopIndex];
            AssertMsg(trip[stopIndex].departureTime >= previousRound()[stop].arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << stop << ", StopIndex: " << stopIndex << ", Time: " << previousRound()[stop].arrivalTime << ", LastDeparture: " << trip[stopIndex].departureTime << ")!");

            StopIndex parentIndex = stopIndex;
            const StopEvent* firstTrip = data.firstTripOfRoute(route);
            while (stopIndex < tripSize - 1) {
                while ((trip > firstTrip) && ((trip - tripSize + stopIndex)->departureTime >= previousRound()[stop].arrivalTime)) {
                    trip -= tripSize;
                    parentIndex = stopIndex;
                }
                stopIndex++;
                stop = stops[stopIndex];
                profiler.countMetric(METRIC_ROUTE_SEGMENTS);
                if (arrivalByRoute(stop, trip[stopIndex].arrivalTime)) {
                    EarliestArrivalLabel& label = currentRound()[stop];
                    label.parent = stops[parentIndex];
                    label.parentDepartureTime = trip[parentIndex].departureTime;
                    label.usesRoute = true;
                    label.routeId = route;
                }
            }
        }
    }

    inline void relaxInitialTransfers() noexcept {
        initialAndFinalTransfers.template addSource<false>(sourceVertex, sourceDepartureTime, sourceVertex);
        initialAndFinalTransfers.initialUpwardSearch();
        initialAndFinalTransfers.downwardSearchToStops();
        for (const StopId stop : data.stops()) {
            const int arrivalTime = initialAndFinalTransfers.getDistance(0, stop);
            if (arrivalByTransfer(stop, arrivalTime)) {
                EarliestArrivalLabel& label = currentRound()[stop];
                label.parent = sourceVertex;
                label.parentDepartureTime = sourceDepartureTime;
                label.usesRoute = false;
                label.transferId = noEdge;
            }
        }
    }

    inline void relaxIntermediateTransfers() noexcept {
        profiler.startPhase();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        for (const StopId stop : stopsUpdatedByRoute) {
            const int earliestArrivalTime = currentRound()[stop].arrivalTime;
            for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
                const StopId toStop = StopId(data.transferGraph.get(ToVertex, edge));
                profiler.countMetric(METRIC_EDGES);
                const int arrivalTime = earliestArrivalTime + data.transferGraph.get(TravelTime, edge);
                AssertMsg(data.isStop(data.transferGraph.get(ToVertex, edge)), "Graph contains edges to non stop vertices!");
                if (arrivalByTransfer(toStop, arrivalTime)) {
                    EarliestArrivalLabel& label = currentRound()[toStop];
                    label.parent = stop;
                    label.parentDepartureTime = earliestArrivalTime;
                    label.usesRoute = false;
                    label.transferId = edge;
                }
            }
            stopsUpdatedByTransfer.insert(stop);
            initialAndFinalTransfers.template addSource<true>(stop, earliestArrivalTime, stop);
        }
        profiler.donePhase(PHASE_TRANSFERS);
        profiler.startPhase();
        initialAndFinalTransfers.relaxFinalTransfers();
        profiler.donePhase(PHASE_FINAL_TRANSFERS);
    }

    inline Round& currentRound() noexcept {
        AssertMsg(!rounds.empty(), "Cannot return current round, because no round exists!");
        return rounds.back();
    }

    inline Round& previousRound() noexcept {
        AssertMsg(rounds.size() >= 2, "Cannot return previous round, because less than two rounds exist!");
        return rounds[rounds.size() - 2];
    }

    inline void startNewRound() noexcept {
        rounds.emplace_back(earliestArrival.size());
        initialAndFinalTransfers.startNewRound();
    }

    inline bool arrivalByRoute(const StopId stop, const int time) noexcept {
        AssertMsg(data.isStop(stop), "Stop " << stop << " is out of range!");
        if (earliestArrival[stop] <= time) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        currentRound()[stop].arrivalTime = time;
        earliestArrival[stop] = time;
        stopsUpdatedByRoute.insert(stop);
        return true;
    }

    inline bool arrivalByTransfer(const StopId stop, const int time) noexcept {
        AssertMsg(data.isStop(stop), "Stop " << stop << " is out of range!");
        if (earliestArrival[stop] <= time) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        currentRound()[stop].arrivalTime = time;
        earliestArrival[stop] = time;
        stopsUpdatedByTransfer.insert(stop);
        return true;
    }

private:
    const Data& data;

    InitialAndFinalTransfers initialAndFinalTransfers;

    Vertex sourceVertex;
    int sourceDepartureTime;
    const IndexedSet<false, Vertex> targetVertices;

    std::vector<Round> rounds;

    std::vector<int> earliestArrival;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Profiler profiler;
};

}
