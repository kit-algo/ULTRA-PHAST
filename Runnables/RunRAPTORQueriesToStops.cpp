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
#include "../DataStructures/RAPTOR/Data.h"
#include "../Helpers/String/String.h"
#include "../Helpers/Timer.h"

using MCR = RAPTOR::OneToAllMCR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler>;

inline void usage() noexcept {
    std::cout << "Usage: RunRAPTORQueriesToStops <MCR binary> <MCR CH data> <UP-RAPTOR binary> <UP-CH data> <number of queries> <seed> <initial transfers: Bucket/RPHAST> <vertex order: DFS/Level> <grouped sweeps: 4/6/8>" << std::endl;
    exit(0);
}

template<bool USE_STOP_BUCKETS, bool USE_DFS_ORDER, size_t GROUPED_SWEEPS>
inline void run(char** argv) noexcept {
    using UPRAPTOR = RAPTOR::UPRAPTOR<GROUPED_SWEEPS, USE_DFS_ORDER, RAPTOR::AggregateProfiler>;

    const std::string mcrFile = argv[1];
    RAPTOR::Data mcrData = RAPTOR::Data::FromBinary(mcrFile);
    mcrData.useImplicitDepartureBufferTimes();
    mcrData.printInfo();
    const std::string mcrCHFile = argv[2];
    CH::CH mcrCH(mcrCHFile);
    const std::string upRAPTORFile = argv[3];
    RAPTOR::Data upRAPTORData = RAPTOR::Data::FromBinary(upRAPTORFile);
    upRAPTORData.useImplicitDepartureBufferTimes();
    upRAPTORData.printInfo();
    const std::string upCHFile = argv[4];
    CH::CH upCH(upCHFile);
    TransferGraph mcrTransferGraph = mcrData.transferGraph;

    const size_t numberOfQueries = String::lexicalCast<size_t>(argv[5]);
    const size_t seed = String::lexicalCast<size_t>(argv[6]);
    srand(seed);

    IndexedSet<false, Vertex> targetSet(upCH.numVertices());
    for (const StopId stop : mcrData.stops()) {
        targetSet.insert(stop);
    }
    Order order = UPRAPTOR::vertexOrder(upCH).splitAt(mcrData.numberOfStops());
    mcrCH.applyVertexOrder(order);
    mcrData.applyVertexOrder(order);
    MCR mcr(mcrData, mcrCH);
    Timer timer;
    UPRAPTOR upRAPTOR = UPRAPTOR::Reordered(upRAPTORData, mcrTransferGraph, upCH, targetSet);
    const double upRAPTORBuildTime = timer.elapsedMicroseconds();

    std::vector<Vertex> sources;
    std::vector<int> departureTimes;

    for (size_t i = 0; i < numberOfQueries; i++) {
        sources.emplace_back(rand() % upCH.numVertices());
        departureTimes.emplace_back(rand() % 24 * 60 * 60);
    }

    for (size_t i = 0; i < sources.size(); i++) {
        mcr.run(sources[i], departureTimes[i]);
    }

    for (size_t i = 0; i < sources.size(); i++) {
        upRAPTOR.run(sources[i], departureTimes[i]);
    }

    std::cout << "UP-RAPTOR setup:" << std::endl;
    std::cout << "\tTime: " << String::musToString(upRAPTORBuildTime) << std::endl;
    std::cout << "\tUpward sweep graph vertices: " << String::prettyInt(upRAPTOR.getUpwardSweepGraphVertices()) << std::endl;
    std::cout << "\tUpward sweep graph edges: " << String::prettyInt(upRAPTOR.getUpwardSweepGraphEdges()) << std::endl;
    std::cout << "\tStop graph vertices: " << String::prettyInt(upRAPTOR.getStopGraphVertices()) << std::endl;
    std::cout << "\tStop graph edges: " << String::prettyInt(upRAPTOR.getStopGraphEdges()) << std::endl;
    std::cout << "\tTarget graph vertices: " << String::prettyInt(upRAPTOR.getTargetGraphVertices()) << std::endl;
    std::cout << "\tTarget graph edges: " << String::prettyInt(upRAPTOR.getTargetGraphEdges()) << std::endl;
    std::cout << std::endl;

    std::cout << "MCR statistics:" << std::endl;
    mcr.getProfiler().printStatistics();
    std::cout << "UP-RAPTOR statistics:" << std::endl;
    upRAPTOR.getProfiler().printStatistics();
}

template<bool USE_STOP_BUCKETS, bool USE_DFS_ORDER>
inline void chooseGroupedSweeps(char** argv) noexcept {
    const size_t groupedSweeps = String::lexicalCast<size_t>(argv[9]);
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
    const std::string orderType = argv[8];
    if (orderType == "DFS") {
        chooseGroupedSweeps<USE_STOP_BUCKETS, true>(argv);
    } else if (orderType == "Level") {
        chooseGroupedSweeps<USE_STOP_BUCKETS, false>(argv);
    } else usage();
}

int main(int argc, char** argv) {
    if (argc < 10) usage();
    const std::string initialTransferType = argv[7];
    if (initialTransferType == "Bucket") {
        chooseOrder<false>(argv);
    } else if (initialTransferType == "RPHAST") {
        chooseOrder<true>(argv);
    } else usage();
    return 0;
}
