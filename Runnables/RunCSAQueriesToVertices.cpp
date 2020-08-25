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
#include "../Helpers/Timer.h"

using MCSA = CSA::OneToAllMCSA<CSA::TransferGraph, CSA::AggregateProfiler>;

inline void usage() noexcept {
    std::cout << "Usage: RunCSAQueriesToVertices <MCSA binary> <UP-CSA binary> <CH data> <number of queries> <seed> <initial transfers: Bucket/RPHAST> <vertex order: DFS/Level>" << std::endl;
    exit(0);
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
    const std::string chFile = argv[3];
    CH::CH ch(chFile);

    const size_t numberOfQueries = String::lexicalCast<size_t>(argv[4]);
    const size_t seed = String::lexicalCast<size_t>(argv[5]);
    srand(seed);

    IndexedSet<false, Vertex> targetSet(Construct::Complete, mCSAData.transferGraph.numVertices());
    mCSAData.applyVertexOrder(UPCSA::vertexOrder(ch));
    MCSA mCSA(mCSAData);
    Timer timer;
    UPCSA upCSA = UPCSA::Reordered(upCSAData, ch, targetSet);
    const double upCSABuildTime = timer.elapsedMicroseconds();

    std::vector<Vertex> sources;
    std::vector<int> departureTimes;

    for (size_t i = 0; i < numberOfQueries; i++) {
        sources.emplace_back(rand() % mCSAData.transferGraph.numVertices());
        departureTimes.emplace_back(rand() % 24 * 60 * 60);
    }

    for (size_t i = 0; i < sources.size(); i++) {
        mCSA.run(sources[i], departureTimes[i]);
    }

    for (size_t i = 0; i < sources.size(); i++) {
        upCSA.run(sources[i], departureTimes[i]);
    }

    std::cout << "UP-CSA setup:" << std::endl;
    std::cout << "\tTime: " << String::musToString(upCSABuildTime) << std::endl;
    std::cout << "\tUpward sweep graph vertices: " << String::prettyInt(upCSA.getUpwardSweepGraphVertices()) << std::endl;
    std::cout << "\tUpward sweep graph edges: " << String::prettyInt(upCSA.getUpwardSweepGraphEdges()) << std::endl;
    std::cout << "\tStop graph vertices: " << String::prettyInt(upCSA.getStopGraphVertices()) << std::endl;
    std::cout << "\tStop graph edges: " << String::prettyInt(upCSA.getStopGraphEdges()) << std::endl;
    std::cout << "\tTarget graph vertices: " << String::prettyInt(upCSA.getTargetGraphVertices()) << std::endl;
    std::cout << "\tTarget graph edges: " << String::prettyInt(upCSA.getTargetGraphEdges()) << std::endl;
    std::cout << std::endl;

    std::cout << "MCSA statistics:" << std::endl;
    mCSA.getProfiler().printStatistics();
    std::cout << "UP-CSA statistics:" << std::endl;
    upCSA.getProfiler().printStatistics();

}

template<bool USE_STOP_BUCKETS>
inline void chooseOrder(char** argv) noexcept {
    const std::string orderType = argv[7];
    if (orderType == "DFS") {
        run<USE_STOP_BUCKETS, true>(argv);
    } else if (orderType == "Level") {
        run<USE_STOP_BUCKETS, false>(argv);
    } else usage();
}

int main(int argc, char** argv) {
    if (argc < 8) usage();
    const std::string initialTransferType = argv[6];
    if (initialTransferType == "Bucket") {
        chooseOrder<false>(argv);
    } else if (initialTransferType == "RPHAST") {
        chooseOrder<true>(argv);
    } else usage();
    return 0;
}
