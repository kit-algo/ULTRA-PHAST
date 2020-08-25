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

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../Algorithms/CH/CH.h"
#include "../Algorithms/CSA/OneToAllMCSA.h"
#include "../Algorithms/CSA/Profiler.h"
#include "../Algorithms/CSA/UPCSA.h"
#include "../DataStructures/Container/Set.h"
#include "../DataStructures/CSA/Data.h"
#include "../Helpers/String/String.h"
#include "../Helpers/Vector/Permutation.h"
#include "../Helpers/Timer.h"

using MCSA = CSA::OneToAllMCSA<CHGraph, CSA::AggregateProfiler>;

inline void usage() noexcept {
    std::cout << "Usage: RunCSAQueriesToBall <MCSA binary> <UP-CSA binary> <number of target sets> <number of targets> <ball size factor> <number of sources> <baseline core degree> <stop factor> <target factor> <seed> <initial transfers: Bucket/RPHAST> <vertex order: DFS/Level>" << std::endl;
    exit(0);
}

inline static std::vector<Vertex> createTargetSet(const TransferGraph& graph, const size_t ballSize, const size_t numTargets) noexcept {
    Dijkstra<TransferGraph> dijkstra(graph);
    std::vector<Vertex> ball;
    do {
        ball.clear();
        const Vertex ballCenter(rand() % graph.numVertices());
        dijkstra.run(ballCenter, noVertex, [&](const Vertex v) {
            ball.emplace_back(v);
            }, [&]() {
                return ball.size() == ballSize;
            });
    } while(ball.size() != ballSize);
    const Permutation permutation(Construct::Random, ballSize);
    std::vector<Vertex> targets;
    for (size_t t = 0; t < numTargets; t++) {
        targets.emplace_back(ball[permutation[t]]);
    }
    return targets;
}

struct CHData {
    CHData(const CH::CH&& ch, const TransferGraph& originalGraph) :
        ch(std::move(ch)),
        numCoreVertices(0) {
        DynamicTransferGraph tempCore;
        tempCore.addVertices(originalGraph.numVertices());
        tempCore[Coordinates] = originalGraph[Coordinates];
        for (const Vertex vertex : tempCore.vertices()) {
            if (ch.isCoreVertex(vertex)) {
                numCoreVertices++;
                for (const Edge edge : ch.forward.edgesFrom(vertex)) {
                    tempCore.addEdge(vertex, ch.forward.get(ToVertex, edge)).set(TravelTime, ch.forward.get(Weight, edge));
                }
            }
        }
        Graph::move(std::move(tempCore), core);
    }

    TransferGraph core;
    CH::CH ch;
    size_t numCoreVertices;
};

inline constexpr int ShortcutWeight = 1024;
inline constexpr int LevelWeight = 256;
inline constexpr int DegreeWeight = 0;

using WitnessSearch = CH::BidirectionalWitnessSearch<CHCoreGraph, 200>;
using GreedyKey = CH::GreedyKey<WitnessSearch>;
using PartialKey = CH::PartialKey<WitnessSearch, GreedyKey>;
using StaggeredKey = CH::StaggeredKey<WitnessSearch, GreedyKey>;
using StopCriterion = CH::CoreCriterion;
using NoStopCriterion = CH::NoStopCriterion;
template<typename KEY_FUNCTION, typename STOP_CRITERION>
using CHBuilder = CH::Builder<WitnessSearch, KEY_FUNCTION, STOP_CRITERION, false, false>;

inline GreedyKey getGreedyKey() noexcept {
    return GreedyKey(ShortcutWeight, LevelWeight, DegreeWeight);
}

template<typename GRAPH, typename KEY_FUNCTION, typename STOP_CRITERION>
inline static CHData buildCH(const GRAPH& originalGraph, const KEY_FUNCTION& keyFunction, const STOP_CRITERION& stopCriterion) noexcept {
    TravelTimeGraph graph;
    Graph::copy(originalGraph, graph);
    Graph::printInfo(graph);
    CHBuilder<KEY_FUNCTION, STOP_CRITERION> chBuilder(std::move(graph), graph[TravelTime], keyFunction, stopCriterion);
    chBuilder.run();
    chBuilder.copyCoreToCH();
    return CHData(CH::CH(std::move(chBuilder)), originalGraph);
}

inline static CHData buildCoreCH(const CSA::Data& csaData, const std::vector<bool>& contractable, const size_t coreDegree) noexcept {
    PartialKey keyFunction(contractable, csaData.transferGraph.numVertices(), getGreedyKey());
    StopCriterion stopCriterion(csaData.numberOfStops(), coreDegree);
    return buildCH(csaData.transferGraph, keyFunction, stopCriterion);
}


inline CH::CH buildULTRACH(const CSA::Data& csaData, const std::vector<Vertex>& targets, const double stopFactor, const double targetFactor) noexcept {
    const size_t stopLimit = csaData.numberOfStops() * stopFactor;
    const size_t targetLimit = targets.size() * targetFactor;
    const size_t targetRound = (targetLimit < stopLimit) ? 2 : 1;
    const size_t stopRound = (targetLimit < stopLimit) ? 1 : 2;
    std::vector<size_t> coreSizes(2);
    coreSizes[targetRound - 1] = targetLimit;
    coreSizes[stopRound - 1] = stopLimit;

    std::vector<size_t> firstContractableRound(csaData.transferGraph.numVertices(), 0);
    for (const Vertex target : targets) {
        firstContractableRound[target] = targetRound;
    }
    for (const StopId stop : csaData.stops()) {
        firstContractableRound[stop] = std::max(firstContractableRound[stop], stopRound);
    }

    StaggeredKey keyFunction(firstContractableRound, coreSizes);
    return buildCH(csaData.transferGraph, keyFunction, NoStopCriterion()).ch;
}


struct BenchmarkData {
    BenchmarkData() :
        baselineCoreBuildTime(0.0),
        baselineCoreVertices(0),
        baselineCoreEdges(0),
        baselineCHEdges(0),
        ultraCHBuildTime(0.0),
        ultraCHEdges(0),
        upCSABuildTime(0.0),
        upwardSweepGraphVertices(0),
        upwardSweepGraphEdges(0),
        stopGraphVertices(0),
        stopGraphEdges(0),
        targetGraphVertices(0),
        targetGraphEdges(0) {
        mCSAProfiler.registerPhases({CSA::PHASE_CLEAR, CSA::PHASE_INITIALIZATION, CSA::PHASE_CONNECTION_SCAN, CSA::PHASE_FINAL_TRANSFERS});
        mCSAProfiler.registerMetrics({CSA::METRIC_CONNECTIONS, CSA::METRIC_EDGES, CSA::METRIC_STOPS_BY_TRIP, CSA::METRIC_STOPS_BY_TRANSFER});
        upCSAProfiler.registerPhases({CSA::PHASE_CLEAR, CSA::PHASE_INITIALIZATION, CSA::PHASE_CONNECTION_SCAN, CSA::PHASE_UPWARD_SWEEP, CSA::PHASE_DOWNWARD_SEARCH});
        upCSAProfiler.registerMetrics({CSA::METRIC_CONNECTIONS, CSA::METRIC_EDGES, CSA::METRIC_STOPS_BY_TRIP, CSA::METRIC_STOPS_BY_TRANSFER});
    }

    inline void print(const size_t numQueries) noexcept {
        std::cout << "Baseline core computation:" << std::endl;
        std::cout << "\tTime: " << String::musToString(baselineCoreBuildTime/numQueries) << std::endl;
        std::cout << "\tCore vertices: " << String::prettyDouble(baselineCoreVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tCore edges: " << String::prettyDouble(baselineCoreEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tCH edges: " << String::prettyDouble(baselineCHEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << std::endl;

        std::cout << "ULTRA CH computation:" << std::endl;
        std::cout << "\tTime: " << String::musToString(ultraCHBuildTime/numQueries) << std::endl;
        std::cout << "\tCH edges: " << String::prettyDouble(ultraCHEdges/(double)numQueries) << std::endl;
        std::cout << std::endl;

        std::cout << "UP-CSA setup:" << std::endl;
        std::cout << "\tTime: " << String::musToString(upCSABuildTime/numQueries) << std::endl;
        std::cout << "\tUpward sweep graph vertices: " << String::prettyDouble(upwardSweepGraphVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tUpward sweep graph edges: " << String::prettyDouble(upwardSweepGraphEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tStop graph vertices: " << String::prettyDouble(stopGraphVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tStop graph edges: " << String::prettyDouble(stopGraphEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tTarget graph vertices: " << String::prettyDouble(targetGraphVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tTarget graph edges: " << String::prettyDouble(targetGraphEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << std::endl;

        std::cout << "MCSA statistics:" << std::endl;
        mCSAProfiler.printStatistics();
        std::cout << "UP-CSA statistics:" << std::endl;
        upCSAProfiler.printStatistics();
    }

    CSA::AggregateProfiler mCSAProfiler;
    CSA::AggregateProfiler upCSAProfiler;
    double baselineCoreBuildTime;
    long long baselineCoreVertices;
    long long baselineCoreEdges;
    long long baselineCHEdges;
    double ultraCHBuildTime;
    long long ultraCHEdges;
    double upCSABuildTime;
    long long upwardSweepGraphVertices;
    long long upwardSweepGraphEdges;
    long long stopGraphVertices;
    long long stopGraphEdges;
    long long targetGraphVertices;
    long long targetGraphEdges;
};

template<typename UPCSA>
inline void runForTargetSet(const CSA::Data& mCSAData, CSA::Data upCSAData, BenchmarkData& data, std::vector<Vertex>& targets, const size_t baselineCoreDegree, const double stopFactor, const double targetFactor, const size_t numSources) noexcept {
    std::vector<bool> isCoreContractable(mCSAData.transferGraph.numVertices(), true);
    for (const Vertex target : targets) {
        isCoreContractable[target] = false;
    }
    for (const StopId stop : mCSAData.stops()) {
        isCoreContractable[stop] = false;
    }

    std::cout << "Building core CH for MCSA" << std::endl;
    Timer timer;
    CHData coreCHData = buildCoreCH(mCSAData, isCoreContractable, baselineCoreDegree);
    data.baselineCoreBuildTime += timer.elapsedMicroseconds();
    data.baselineCoreVertices += coreCHData.numCoreVertices;
    data.baselineCoreEdges += coreCHData.core.numEdges();
    data.baselineCHEdges += coreCHData.ch.numEdges();
    CSA::Data coreCSAData = mCSAData;
    coreCSAData.transferGraph = coreCHData.core;
    CHGraph upwardGraph = coreCHData.ch.forward;
    std::cout << "Building CH for UP-CSA" << std::endl;
    timer.restart();
    CH::CH ultraCH = buildULTRACH(mCSAData, targets, stopFactor, targetFactor);
    data.ultraCHBuildTime += timer.elapsedMicroseconds();
    data.ultraCHEdges += ultraCH.numEdges();
    IndexedSet<false, Vertex> targetSet(mCSAData.transferGraph.numVertices(), targets);

    std::cout << "Reordering vertices" << std::endl;
    const Order order = UPCSA::vertexOrder(ultraCH).splitAt(mCSAData.numberOfStops());
    coreCSAData.applyVertexOrder(order);
    upwardGraph.applyVertexOrder(order);

    std::cout << "Building algorithms" << std::endl;
    MCSA mCSA(coreCSAData, upwardGraph, Weight);
    timer.restart();
    UPCSA upCSA = UPCSA::Reordered(upCSAData, ultraCH, targetSet);
    data.upCSABuildTime += timer.elapsedMicroseconds();
    data.upwardSweepGraphVertices += upCSA.getUpwardSweepGraphVertices();
    data.upwardSweepGraphEdges += upCSA.getUpwardSweepGraphEdges();
    data.stopGraphVertices += upCSA.getStopGraphVertices();
    data.stopGraphEdges += upCSA.getStopGraphEdges();
    data.targetGraphVertices += upCSA.getTargetGraphVertices();
    data.targetGraphEdges += upCSA.getTargetGraphEdges();

    std::vector<Vertex> sources;
    std::vector<int> departureTimes;

    for (size_t j = 0; j < numSources; j++) {
        sources.emplace_back(rand() % ultraCH.numVertices());
        departureTimes.emplace_back(rand() % 24 * 60 * 60);
    }

    std::cout << "Running queries" << std::endl;
    for (size_t j = 0; j < sources.size(); j++) {
        mCSA.run(sources[j], departureTimes[j]);
    }

    for (size_t j = 0; j < sources.size(); j++) {
        upCSA.run(sources[j], departureTimes[j]);
    }

    data.mCSAProfiler += mCSA.getProfiler();
    data.upCSAProfiler += upCSA.getProfiler();
}

template<bool USE_STOP_BUCKETS, bool USE_DFS_ORDER>
inline void run(char** argv) noexcept {
    using UPCSA = CSA::UPCSA<USE_STOP_BUCKETS, USE_DFS_ORDER, CSA::AggregateProfiler>;

    const std::string mCSAFile = argv[1];
    CSA::Data mCSAData = CSA::Data::FromBinary(mCSAFile);
    mCSAData.sortConnectionsAscending();
    mCSAData.printInfo();
    const std::string upCSAFile = argv[2];
    CSA::Data upCSAData = CSA::Data::FromBinary(upCSAFile);
    upCSAData.sortConnectionsAscending();
    upCSAData.printInfo();

    const size_t numTargetSets = String::lexicalCast<size_t>(argv[3]);
    const size_t numTargets = String::lexicalCast<size_t>(argv[4]);
    const double ballSizeFactor = String::lexicalCast<double>(argv[5]);
    const size_t ballSize = std::min(mCSAData.transferGraph.numVertices(), size_t(numTargets * ballSizeFactor));
    const size_t numSources = String::lexicalCast<size_t>(argv[6]);
    const size_t baselineCoreDegree = String::lexicalCast<size_t>(argv[7]);
    const double stopFactor = String::lexicalCast<double>(argv[8]);
    const double targetFactor = String::lexicalCast<double>(argv[9]);
    const size_t seed = String::lexicalCast<size_t>(argv[10]);
    srand(seed);

    BenchmarkData data;
    for (size_t i = 0; i < numTargetSets; i++) {
        std::cout << "Creating target set " << i << std::endl;
        std::vector<Vertex> targets = createTargetSet(mCSAData.transferGraph, ballSize, numTargets);
        runForTargetSet<UPCSA>(mCSAData, upCSAData, data, targets, baselineCoreDegree, stopFactor, targetFactor, numSources);
    }
    data.print(numTargetSets);
}

template<bool USE_STOP_BUCKETS>
inline void chooseOrder(char** argv) noexcept {
    const std::string orderType = argv[12];
    if (orderType == "DFS") {
        run<USE_STOP_BUCKETS, true>(argv);
    } else if (orderType == "Level") {
        run<USE_STOP_BUCKETS, false>(argv);
    } else usage();
}

int main(int argc, char** argv) {
    if (argc < 13) usage();
    const std::string initialTransferType = argv[11];
    if (initialTransferType == "Bucket") {
        chooseOrder<false>(argv);
    } else if (initialTransferType == "RPHAST") {
        chooseOrder<true>(argv);
    } else usage();
    return 0;
}
