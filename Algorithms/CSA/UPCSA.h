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

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "../CH/CH.h"
#include "../CH/CHUtils.h"
#include "../CH/UPQuery/UPQuery.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "Profiler.h"

namespace CSA {

template<bool USE_STOP_BUCKETS, bool USE_DFS_ORDER, typename PROFILER = NoProfiler>
class UPCSA {

public:
    constexpr static bool UseStopBuckets = USE_STOP_BUCKETS;
    constexpr static bool UseDFSOrder = USE_DFS_ORDER;
    using Profiler = PROFILER;
    constexpr static bool Debug = Meta::Equals<Profiler, SimpleProfiler>();
    using Type = UPCSA<UseStopBuckets, UseDFSOrder, Profiler>;
    using InitialAndFinalTransfers = CH::UPQuery<UseStopBuckets, true, Debug>;
    using TripFlag = ConnectionId;

private:
    struct ParentLabel {
        ParentLabel(const Vertex parent = noVertex, const bool reachedByTransfer = false, const TripId tripId = noTripId) :
            parent(parent),
            reachedByTransfer(reachedByTransfer),
            tripId(tripId) {
        }

        Vertex parent;
        bool reachedByTransfer;
        union {
            TripId tripId;
            Edge transferId;
        };
    };

    UPCSA(const Data& data, const CH::CH& chData, const Order&& chOrder, const IndexedSet<false, Vertex>& targets, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialAndFinalTransfers(chData, std::move(chOrder), data.numberOfStops(), targets),
        sourceVertex(noVertex),
        sourceDepartureTime(never),
        targetVertices(targets),
        tripReached(data.numberOfTrips(), TripFlag()),
        arrivalTime(data.numberOfStops(), never),
        arrivalTimeByTrip(data.numberOfStops(), never),
        parentLabel(data.numberOfStops()),
        profiler(profilerTemplate) {
        AssertMsg(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        AssertMsg(!Graph::hasLoops(data.transferGraph), "Shortcut graph may not have loops!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN, PHASE_UPWARD_SWEEP, PHASE_DOWNWARD_SEARCH});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
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

    inline static UPCSA Reordered(Data& data, CH::CH& chData, IndexedSet<false, Vertex>& targets, const Profiler& profilerTemplate = Profiler()) noexcept {
        Order chOrder = vertexOrder(chData);
        Order fullOrder = chOrder.splitAt(data.numberOfStops());
        Order stopOrder = fullOrder;
        stopOrder.resize(data.numberOfStops());

        data.applyStopOrder(stopOrder);
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

        return UPCSA(data, chData, std::move(phastOrder), targets, profilerTemplate);
    }

    inline void run(const Vertex source, const int departureTime) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceVertex = source;
        sourceDepartureTime = departureTime;
        initialAndFinalTransfers.initialize();
        if (data.isStop(source)) {
            arrivalTime[source] = departureTime;
            arrivalTimeByTrip[source] = departureTime;
        }
        runInitialTransfers();
        const ConnectionId firstConnection = firstReachableConnection(departureTime);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections(firstConnection, ConnectionId(data.connections.size()));
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.startPhase();
        for (const StopId stop : data.stops()) {
            initialAndFinalTransfers.template addSource<true>(stop, arrivalTime[stop], stop);
        }
        initialAndFinalTransfers.upwardSweep();
        profiler.donePhase(PHASE_UPWARD_SWEEP);
        profiler.startPhase();
        initialAndFinalTransfers.downwardSearchToTargets();
        profiler.donePhase(PHASE_DOWNWARD_SEARCH);

        profiler.done();
    }

    inline bool reachable(const Vertex vertex) noexcept {
        AssertMsg(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        return getEarliestArrivalTime(vertex) < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) noexcept {
        AssertMsg(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        return initialAndFinalTransfers.getDistance(vertex);
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
    inline void clear() {
        sourceVertex = noVertex;
        sourceDepartureTime = never;
        Vector::fill(arrivalTime, never);
        Vector::fill(arrivalTimeByTrip, never);
        Vector::fill(tripReached, TripFlag());
        Vector::fill(parentLabel, ParentLabel());
    }

    inline ConnectionId firstReachableConnection(const int departureTime) const noexcept {
        return ConnectionId(Vector::lowerBound(data.connections, departureTime, [](const Connection& connection, const int time) {
            return connection.departureTime < time;
        }));
    }

    inline void scanConnections(const ConnectionId begin, const ConnectionId end) noexcept {
        for (ConnectionId i = begin; i < end; i++) {
            const Connection& connection = data.connections[i];
            if (connectionIsReachable(connection, i)) {
                profiler.countMetric(METRIC_CONNECTIONS);
                arrivalByTrip(connection.arrivalStopId, connection.arrivalTime, connection.tripId);
            }
        }
    }

    inline bool connectionIsReachableFromStop(const Connection& connection) const noexcept {
        return arrivalTime[connection.departureStopId] <= connection.departureTime - data.minTransferTime(connection.departureStopId);
    }

    inline bool connectionIsReachableFromTrip(const Connection& connection) const noexcept {
        return tripReached[connection.tripId] != TripFlag();
    }

    inline bool connectionIsReachable(const Connection& connection, const ConnectionId id) noexcept {
        if (connectionIsReachableFromTrip(connection)) return true;
        if (connectionIsReachableFromStop(connection)) {
            tripReached[connection.tripId] = id;
            return true;
        }
        return false;
    }

    inline void arrivalByTrip(const StopId stop, const int time, const TripId trip) noexcept {
        if (arrivalTimeByTrip[stop] <= time) return;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        arrivalTimeByTrip[stop] = time;

        for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
            profiler.countMetric(METRIC_EDGES);
            const StopId toStop = StopId(data.transferGraph.get(ToVertex, edge));
            const int newArrivalTime = time + data.transferGraph.get(TravelTime, edge);
            arrivalByTransfer(toStop, newArrivalTime, stop, edge);
        }

        if (arrivalTime[stop] <= time) return;
        arrivalTime[stop] = time;
        parentLabel[stop].parent = data.connections[tripReached[trip]].departureStopId;
        parentLabel[stop].reachedByTransfer = false;
        parentLabel[stop].tripId = trip;
    }

    inline void runInitialTransfers() noexcept {
        initialAndFinalTransfers.template addSource<false>(sourceVertex, sourceDepartureTime, sourceVertex);
        initialAndFinalTransfers.upwardSearch();
        initialAndFinalTransfers.downwardSearchToStops();
        for (const StopId stop : data.stops()) {
            arrivalByTransfer(stop, initialAndFinalTransfers.getDistance(stop), sourceVertex, noEdge);
        }
    }

    inline void arrivalByTransfer(const StopId stop, const int time, const Vertex parent, const Edge edge) noexcept {
        if (arrivalTime[stop] <= time) return;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        arrivalTime[stop] = time;
        parentLabel[stop].parent = parent;
        parentLabel[stop].reachedByTransfer = true;
        parentLabel[stop].transferId = edge;
    }

private:
    const Data& data;

    InitialAndFinalTransfers initialAndFinalTransfers;

    Vertex sourceVertex;
    int sourceDepartureTime;
    const IndexedSet<false, Vertex> targetVertices;

    std::vector<TripFlag> tripReached;
    std::vector<int> arrivalTime;
    std::vector<int> arrivalTimeByTrip;
    std::vector<ParentLabel> parentLabel;

    Profiler profiler;

};
}
