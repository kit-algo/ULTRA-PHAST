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
#include "../Algorithms/RAPTOR/OneToAllMCR.h"
#include "../Algorithms/RAPTOR/Profiler.h"
#include "../Algorithms/RAPTOR/UPRAPTOR.h"
#include "../DataStructures/Container/Set.h"
#include "../DataStructures/CSA/Data.h"
#include "../Helpers/String/String.h"
#include "../Helpers/Vector/Permutation.h"
#include "../Helpers/Timer.h"

using MCR = RAPTOR::OneToAllMCR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler>;

inline void usage() noexcept {
    std::cout << "Usage: RunRAPTORQueriesToBall <MCR binary> <UP-RAPTOR binary> <number of target sets> <number of targets> <ball size factor> <number of sources> <baseline core degree> <stop factor> <target factor> <seed> <initial transfers: Bucket/RPHAST> <vertex order: DFS/Level> <grouped sweeps: 4/6/8>" << std::endl;
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

inline static CHData buildCoreCH(const RAPTOR::Data& raptorData, const std::vector<bool>& contractable, const size_t coreDegree) noexcept {
    PartialKey keyFunction(contractable, raptorData.transferGraph.numVertices(), getGreedyKey());
    StopCriterion stopCriterion(raptorData.numberOfStops(), coreDegree);
    return buildCH(raptorData.transferGraph, keyFunction, stopCriterion);
}


inline CH::CH buildULTRACH(const RAPTOR::Data& raptorData, const std::vector<Vertex>& targets, const double stopFactor, const double targetFactor) noexcept {
    const size_t stopLimit = raptorData.numberOfStops() * stopFactor;
    const size_t targetLimit = targets.size() * targetFactor;
    const size_t targetRound = (targetLimit < stopLimit) ? 2 : 1;
    const size_t stopRound = (targetLimit < stopLimit) ? 1 : 2;
    std::vector<size_t> coreSizes(2);
    coreSizes[targetRound - 1] = targetLimit;
    coreSizes[stopRound - 1] = stopLimit;

    std::vector<size_t> firstContractableRound(raptorData.transferGraph.numVertices(), 0);
    for (const Vertex target : targets) {
        firstContractableRound[target] = targetRound;
    }
    for (const StopId stop : raptorData.stops()) {
        firstContractableRound[stop] = std::max(firstContractableRound[stop], stopRound);
    }

    StaggeredKey keyFunction(firstContractableRound, coreSizes);
    return buildCH(raptorData.transferGraph, keyFunction, NoStopCriterion()).ch;
}


struct BenchmarkData {
    BenchmarkData() :
        baselineCoreBuildTime(0.0),
        baselineCoreVertices(0),
        baselineCoreEdges(0),
        baselineCHEdges(0),
        ultraCHBuildTime(0.0),
        ultraCHEdges(0),
        upRAPTORBuildTime(0.0),
        upwardSweepGraphVertices(0),
        upwardSweepGraphEdges(0),
        stopGraphVertices(0),
        stopGraphEdges(0),
        targetGraphVertices(0),
        targetGraphEdges(0) {
        mcrProfiler.registerExtraRounds({RAPTOR::EXTRA_ROUND_CLEAR, RAPTOR::EXTRA_ROUND_INITIALIZATION});
        mcrProfiler.registerPhases({RAPTOR::PHASE_INITIALIZATION, RAPTOR::PHASE_COLLECT, RAPTOR::PHASE_SCAN, RAPTOR::PHASE_TRANSFERS});
        mcrProfiler.registerMetrics({RAPTOR::METRIC_ROUTES, RAPTOR::METRIC_ROUTE_SEGMENTS, RAPTOR::METRIC_VERTICES, RAPTOR::METRIC_EDGES, RAPTOR::METRIC_STOPS_BY_TRIP, RAPTOR::METRIC_STOPS_BY_TRANSFER});
        upRAPTORProfiler.registerExtraRounds({RAPTOR::EXTRA_ROUND_CLEAR, RAPTOR::EXTRA_ROUND_INITIALIZATION, RAPTOR::EXTRA_ROUND_FINAL_TRANSFERS});
        upRAPTORProfiler.registerPhases({RAPTOR::PHASE_INITIALIZATION, RAPTOR::PHASE_COLLECT, RAPTOR::PHASE_SCAN, RAPTOR::PHASE_TRANSFERS, RAPTOR::PHASE_FINAL_TRANSFERS});
        upRAPTORProfiler.registerMetrics({RAPTOR::METRIC_ROUTES, RAPTOR::METRIC_ROUTE_SEGMENTS, RAPTOR::METRIC_EDGES, RAPTOR::METRIC_STOPS_BY_TRIP, RAPTOR::METRIC_STOPS_BY_TRANSFER});
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

        std::cout << "UP-RAPTOR setup:" << std::endl;
        std::cout << "\tTime: " << String::musToString(upRAPTORBuildTime/numQueries) << std::endl;
        std::cout << "\tUpward sweep graph vertices: " << String::prettyDouble(upwardSweepGraphVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tUpward sweep graph edges: " << String::prettyDouble(upwardSweepGraphEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tStop graph vertices: " << String::prettyDouble(stopGraphVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tStop graph edges: " << String::prettyDouble(stopGraphEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tTarget graph vertices: " << String::prettyDouble(targetGraphVertices/static_cast<double>(numQueries)) << std::endl;
        std::cout << "\tTarget graph edges: " << String::prettyDouble(targetGraphEdges/static_cast<double>(numQueries)) << std::endl;
        std::cout << std::endl;

        std::cout << "MCR statistics:" << std::endl;
        mcrProfiler.printStatistics();
        std::cout << "UP-RAPTOR statistics:" << std::endl;
        upRAPTORProfiler.printStatistics();
    }

    RAPTOR::AggregateProfiler mcrProfiler;
    RAPTOR::AggregateProfiler upRAPTORProfiler;
    double baselineCoreBuildTime;
    long long baselineCoreVertices;
    long long baselineCoreEdges;
    long long baselineCHEdges;
    double ultraCHBuildTime;
    long long ultraCHEdges;
    double upRAPTORBuildTime;
    long long upwardSweepGraphVertices;
    long long upwardSweepGraphEdges;
    long long stopGraphVertices;
    long long stopGraphEdges;
    long long targetGraphVertices;
    long long targetGraphEdges;
};

template<typename UPRAPTOR>
inline void runForTargetSet(const RAPTOR::Data& mcrData, RAPTOR::Data upRAPTORData, BenchmarkData& data, std::vector<Vertex>& targets, const size_t baselineCoreDegree, const double stopFactor, const double targetFactor, const size_t numSources) noexcept {
    std::vector<bool> isCoreContractable(mcrData.transferGraph.numVertices(), true);
    for (const Vertex target : targets) {
        isCoreContractable[target] = false;
    }
    for (const StopId stop : mcrData.stops()) {
        isCoreContractable[stop] = false;
    }

    std::cout << "Building core CH for MCR" << std::endl;
    Timer timer;
    CHData coreCHData = buildCoreCH(mcrData, isCoreContractable, baselineCoreDegree);
    data.baselineCoreBuildTime += timer.elapsedMicroseconds();
    data.baselineCoreVertices += coreCHData.numCoreVertices;
    data.baselineCoreEdges += coreCHData.core.numEdges();
    data.baselineCHEdges += coreCHData.ch.numEdges();
    RAPTOR::Data coreRAPTORData = mcrData;
    coreRAPTORData.transferGraph = coreCHData.core;
    CH::CH coreCH = coreCHData.ch;
    std::cout << "Building CH for UP-RAPTOR" << std::endl;
    timer.restart();
    CH::CH ultraCH = buildULTRACH(mcrData, targets, stopFactor, targetFactor);
    data.ultraCHBuildTime += timer.elapsedMicroseconds();
    data.ultraCHEdges += ultraCH.numEdges();
    IndexedSet<false, Vertex> targetSet(mcrData.transferGraph.numVertices(), targets);

    TransferGraph coreGraph = coreRAPTORData.transferGraph;
    std::cout << "Reordering vertices" << std::endl;
    const Order order = UPRAPTOR::vertexOrder(ultraCH).splitAt(mcrData.numberOfStops());
    coreRAPTORData.applyVertexOrder(order);
    coreCH.applyVertexOrder(order);

    std::cout << "Building algorithms" << std::endl;
    MCR mcr(coreRAPTORData, coreCH);
    timer.restart();
    UPRAPTOR upRAPTOR = UPRAPTOR::Reordered(upRAPTORData, coreGraph, ultraCH, targetSet);
    data.upRAPTORBuildTime += timer.elapsedMicroseconds();
    data.upwardSweepGraphVertices += upRAPTOR.getUpwardSweepGraphVertices();
    data.upwardSweepGraphEdges += upRAPTOR.getUpwardSweepGraphEdges();
    data.stopGraphVertices += upRAPTOR.getStopGraphVertices();
    data.stopGraphEdges += upRAPTOR.getStopGraphEdges();
    data.targetGraphVertices += upRAPTOR.getTargetGraphVertices();
    data.targetGraphEdges += upRAPTOR.getTargetGraphEdges();

    std::vector<Vertex> sources;
    std::vector<int> departureTimes;

    for (size_t j = 0; j < numSources; j++) {
        sources.emplace_back(rand() % ultraCH.numVertices());
        departureTimes.emplace_back(rand() % 24 * 60 * 60);
    }

    std::cout << "Running queries" << std::endl;
    for (size_t j = 0; j < sources.size(); j++) {
        mcr.run(sources[j], departureTimes[j]);
    }

    for (size_t j = 0; j < sources.size(); j++) {
        upRAPTOR.run(sources[j], departureTimes[j]);
    }

    data.mcrProfiler += mcr.getProfiler();
    data.upRAPTORProfiler += upRAPTOR.getProfiler();
}

template<bool USE_STOP_BUCKETS, bool USE_DFS_ORDER, size_t GROUPED_SWEEPS>
inline void run(char** argv) noexcept {
    using UPRAPTOR = RAPTOR::UPRAPTOR<GROUPED_SWEEPS, USE_DFS_ORDER, RAPTOR::AggregateProfiler>;

    const std::string mcrFile = argv[1];
    RAPTOR::Data mcrData = RAPTOR::Data::FromBinary(mcrFile);
    mcrData.useImplicitDepartureBufferTimes();
    mcrData.printInfo();
    const std::string upRAPTORFile = argv[2];
    RAPTOR::Data upRAPTORData = RAPTOR::Data::FromBinary(upRAPTORFile);
    upRAPTORData.useImplicitDepartureBufferTimes();
    upRAPTORData.printInfo();

    const size_t numTargetSets = String::lexicalCast<size_t>(argv[3]);
    const size_t numTargets = String::lexicalCast<size_t>(argv[4]);
    const double ballSizeFactor = String::lexicalCast<double>(argv[5]);
    const size_t ballSize = std::min(mcrData.transferGraph.numVertices(), size_t(numTargets * ballSizeFactor));
    const size_t numSources = String::lexicalCast<size_t>(argv[6]);
    const size_t baselineCoreDegree = String::lexicalCast<size_t>(argv[7]);
    const double stopFactor = String::lexicalCast<double>(argv[8]);
    const double targetFactor = String::lexicalCast<double>(argv[9]);
    const size_t seed = String::lexicalCast<size_t>(argv[10]);
    srand(seed);

    BenchmarkData data;
    for (size_t i = 0; i < numTargetSets; i++) {
        std::cout << "Creating target set " << i << std::endl;
        std::vector<Vertex> targets = createTargetSet(mcrData.transferGraph, ballSize, numTargets);
        runForTargetSet<UPRAPTOR>(mcrData, upRAPTORData, data, targets, baselineCoreDegree, stopFactor, targetFactor, numSources);
    }
    data.print(numTargetSets);
}

template<bool USE_STOP_BUCKETS, bool USE_DFS_ORDER>
inline void chooseGroupedSweeps(char** argv) noexcept {
    const size_t groupedSweeps = String::lexicalCast<size_t>(argv[13]);
    if (groupedSweeps == 4) {
        run<USE_STOP_BUCKETS, USE_DFS_ORDER, 4>(argv);
    } else if (groupedSweeps == 6) {
        run<USE_STOP_BUCKETS, USE_DFS_ORDER, 6>(argv);
    } else if (groupedSweeps == 8) {
        run<USE_STOP_BUCKETS, USE_DFS_ORDER, 8>(argv);
    } else usage();
}

template<bool USE_STOP_BUCKETS>
inline void chooseOrder(char** argv) noexcept {
    const std::string orderType = argv[12];
    if (orderType == "DFS") {
        chooseGroupedSweeps<USE_STOP_BUCKETS, true>(argv);
    } else if (orderType == "Level") {
        chooseGroupedSweeps<USE_STOP_BUCKETS, false>(argv);
    } else usage();
}

int main(int argc, char** argv) {
    if (argc < 14) usage();
    const std::string initialTransferType = argv[11];
    if (initialTransferType == "Bucket") {
        chooseOrder<false>(argv);
    } else if (initialTransferType == "RPHAST") {
        chooseOrder<true>(argv);
    } else usage();
    return 0;
}
