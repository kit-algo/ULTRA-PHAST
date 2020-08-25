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

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Container/ExternalKHeap.h"
#include "Profiler.h"

namespace CSA {

template<typename INITIAL_TRANSFER_GRAPH, typename PROFILER = NoProfiler>
class OneToAllMCSA {

public:
    using InitialTransferGraph = INITIAL_TRANSFER_GRAPH;
    using Profiler = PROFILER;
    using Type = OneToAllMCSA<InitialTransferGraph, Profiler>;
    using TripFlag = ConnectionId;

private:
    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel(int* const arrivalTime) :
            ExternalKHeapElement(),
            arrivalTime(arrivalTime) {
        }

        inline int getArrivalTime() const noexcept {
            return *arrivalTime;
        }

        inline bool hasSmallerKey(const DijkstraLabel* other) const noexcept {return getArrivalTime() < other->getArrivalTime();}

        int* const arrivalTime;
    };

public:
    template<typename ATTRIBUTE>
    OneToAllMCSA(const Data& data, const InitialTransferGraph& initialTransferGraph, const ATTRIBUTE weight, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialTransferGraph(initialTransferGraph),
        initialTransferWeight(initialTransferGraph[weight]),
        sourceVertex(noVertex),
        tripReached(data.numberOfTrips(), TripFlag()),
        arrivalTime(initialTransferGraph.numVertices(), never),
        parent(data.numberOfStops(), noVertex),
        parentTrip(data.numberOfStops(), noTripId),
        dijkstraParent(initialTransferGraph.numVertices(), noVertex),
        profiler(profilerTemplate) {
        AssertMsg(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN, PHASE_FINAL_TRANSFERS});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
        for (const Vertex vertex : initialTransferGraph.vertices()) {
            dijkstraLabels.emplace_back(&arrivalTime[vertex]);
        }
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    OneToAllMCSA(const Data& data, const Profiler& profilerTemplate = Profiler()) :
        OneToAllMCSA(data, data.transferGraph, TravelTime, profilerTemplate) {
    }

    inline void run(const Vertex source, const int departureTime) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceVertex = source;
        dijkstraParent[source] = source;
        arrivalByTransfer(source, departureTime, source);
        runInitialTransfers();
        const ConnectionId firstConnection = firstReachableConnection(departureTime);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections(firstConnection, ConnectionId(data.connections.size()));
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.startPhase();
        runDijkstra(INFTY);
        profiler.donePhase(PHASE_FINAL_TRANSFERS);

        profiler.done();
    }

    inline bool reachable(const Vertex vertex) const noexcept {
        return arrivalTime[vertex] < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) const noexcept {
        return arrivalTime[vertex];
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() {
        sourceVertex = noVertex;
        Vector::fill(tripReached, TripFlag());
        Vector::fill(arrivalTime, never);
        Vector::fill(parent, noVertex);
        Vector::fill(parentTrip, noTripId);
        Vector::fill(dijkstraParent, noVertex);
        queue.clear();
    }

    inline void runInitialTransfers() noexcept {
        while(!queue.empty()) {
            settle(initialTransferGraph, initialTransferWeight);
        }
    }

    inline ConnectionId firstReachableConnection(const int departureTime) const noexcept {
        return ConnectionId(Vector::lowerBound(data.connections, departureTime, [](const Connection& connection, const int time) {
            return connection.departureTime < time;
        }));
    }

    inline void scanConnections(const ConnectionId begin, const ConnectionId end) noexcept {
        for (ConnectionId i = begin; i < end; i++) {
            const Connection& connection = data.connections[i];
            runDijkstra(connection.departureTime);
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
        if (arrivalTime[stop] <= time) return;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        arrivalTime[stop] = time;
        queue.update(&dijkstraLabels[stop]);
        parent[stop] = data.connections[tripReached[trip]].departureStopId;
        parentTrip[stop] = trip;
        dijkstraParent[stop] = stop;
    }

    inline void runDijkstra(const int nextDepartureTime) noexcept {
        while ((!queue.empty()) && (queue.min().getArrivalTime() <= nextDepartureTime)) {
            settle(data.transferGraph, data.transferGraph[TravelTime]);
        }
    }

    template<typename GRAPH>
    inline void settle(const GRAPH& graph, const std::vector<int>& weight) noexcept {
        const Vertex u = Vertex(queue.extractFront() - &(dijkstraLabels[0]));
        const int time = arrivalTime[u];
        for (const Edge edge : graph.edgesFrom(u)) {
            profiler.countMetric(METRIC_EDGES);
            const Vertex v = graph.get(ToVertex, edge);
            const int newArrivalTime = time + weight[edge];
            if (newArrivalTime < arrivalTime[v]) {
                arrivalByTransfer(v, newArrivalTime, u);
            }
        }
        if (data.isStop(u)) {
            profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        }
    }

    inline void arrivalByTransfer(const Vertex vertex, const int time, const Vertex parentVertex) noexcept {
        arrivalTime[vertex] = time;
        queue.update(&dijkstraLabels[vertex]);
        dijkstraParent[vertex] = dijkstraParent[parentVertex];
        if (data.isStop(vertex)) {
            parent[vertex] = dijkstraParent[parentVertex];
            parentTrip[vertex] = noTripId;
        }
    }

private:
    const Data& data;
    const InitialTransferGraph& initialTransferGraph;
    const std::vector<int>& initialTransferWeight;

    Vertex sourceVertex;

    std::vector<TripFlag> tripReached;
    std::vector<int> arrivalTime;
    std::vector<Vertex> parent;
    std::vector<TripId> parentTrip;
    std::vector<Vertex> dijkstraParent;

    ExternalKHeap<2, DijkstraLabel> queue;
    std::vector<DijkstraLabel> dijkstraLabels;

    Profiler profiler;

};
}
