/**********************************************************************************

 Copyright (c) 2019 Tobias Zündorf

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

#include <vector>

#include "CH.h"
#include "Preprocessing/CHData.h"

#include "../DepthFirstSearch.h"

#include "../../Helpers/Helpers.h"
#include "../../Helpers/Vector/Vector.h"

namespace CH {

inline std::vector<Vertex> getOrder(const CH& ch) noexcept {
    std::vector<Vertex> order;
    order.reserve(ch.numVertices());
    depthFirstSearch(ch, NoOperation, NoOperation, [&](const Edge edge, const Vertex) {
        order.push_back(ch.get(ToVertex, edge));
    }, NoOperation, [&](const Vertex root) {
        order.push_back(root);
    });
    size_t front = 0;
    size_t back = order.size() - 1;
    std::vector<Vertex> coreOrder(order.size());
    while (!order.empty()) {
        if (ch.isCoreVertex(order.back())) {
            coreOrder[back] = order.back();
            back--;
        } else {
            coreOrder[front] = order.back();
            front++;
        }
        order.pop_back();
    }
    return coreOrder;
}

inline std::vector<uint16_t> vertexLevelTopDown(const CHGraph& ch, const std::vector<Vertex>& order) noexcept {
    std::vector<uint16_t> level(order.size(), 0);
    for (Vertex vertex : descending(order)) {
        for (Edge edge : ch.edgesFrom(vertex)) {
            level[vertex] = std::max<uint16_t>(level[vertex], level[ch.get(ToVertex, edge)] + 1);
        }
    }
    return level;
}

inline std::vector<uint16_t> vertexLevelTopDown(const CH& ch, const std::vector<Vertex>& order) noexcept {
    std::vector<uint16_t> level(order.size(), 0);
    for (Vertex vertex : descending(order)) {
        if (ch.isCoreVertex(vertex)) continue;
        for (Edge edge : ch.forward.edgesFrom(vertex)) {
            level[vertex] = std::max<uint16_t>(level[vertex], level[ch.forward.get(ToVertex, edge)] + 1);
        }
        for (Edge edge : ch.backward.edgesFrom(vertex)) {
            level[vertex] = std::max<uint16_t>(level[vertex], level[ch.backward.get(ToVertex, edge)] + 1);
        }
    }
    return level;
}

inline std::vector<uint16_t> vertexLevelBottomUp(const CH& ch, const std::vector<Vertex>& order) noexcept {
    std::vector<uint16_t> level(order.size(), 0);
    for (Vertex vertex : order) {
        if (ch.isCoreVertex(vertex)) continue;
        for (Edge edge : ch.forward.edgesFrom(vertex)) {
            level[ch.forward.get(ToVertex, edge)] = std::max<uint16_t>(level[vertex] + 1, level[ch.forward.get(ToVertex, edge)]);
        }
        for (Edge edge : ch.backward.edgesFrom(vertex)) {
            level[ch.backward.get(ToVertex, edge)] = std::max<uint16_t>(level[vertex] + 1, level[ch.backward.get(ToVertex, edge)]);
        }
    }
    return level;
}

inline std::vector<std::vector<Vertex>> verticesByLevel(const CH& ch, const std::vector<uint16_t>& level) noexcept {
    std::vector<std::vector<Vertex>> verticesByLevel;
    for (Vertex vertex : ch.vertices()) {
        if (verticesByLevel.size() <= level[vertex]) verticesByLevel.resize(level[vertex] + 1);
        verticesByLevel[level[vertex]].push_back(vertex);
    }
    return verticesByLevel;
}

inline std::vector<size_t> verticesPerLevel(const CH& ch, const std::vector<uint16_t>& level) noexcept {
    std::vector<size_t> verticesPerLevel;
    for (Vertex vertex : ch.vertices()) {
        if (verticesPerLevel.size() <= level[vertex]) verticesPerLevel.resize(level[vertex] + 1, 0);
        verticesPerLevel[level[vertex]]++;
    }
    return verticesPerLevel;
}

inline std::vector<std::vector<Vertex>> verticesByLevelTopDown(const CH& ch, const std::vector<Vertex>& order) noexcept {
    return verticesByLevel(ch, vertexLevelTopDown(ch, order));
}

inline std::vector<size_t> verticesPerLevelTopDown(const CH& ch, const std::vector<Vertex>& order) noexcept {
    return verticesPerLevel(ch, vertexLevelTopDown(ch, order));
}

inline std::vector<Vertex> getLevelOrderTopDown(const CH& ch) noexcept {
    return Vector::flatten(verticesByLevel(ch, vertexLevelTopDown(ch, getOrder(ch))));
}

inline std::vector<std::vector<Vertex>> verticesByLevelBottomUp(const CH& ch, const std::vector<Vertex>& order) noexcept {
    return verticesByLevel(ch, vertexLevelTopDown(ch, order));
}

inline std::vector<size_t> verticesPerLevelBottomUp(const CH& ch, const std::vector<Vertex>& order) noexcept {
    return verticesPerLevel(ch, vertexLevelBottomUp(ch, order));
}

inline std::vector<Vertex> getLevelOrderBottomUp(const CH& ch) noexcept {
    return Vector::flatten(verticesByLevel(ch, vertexLevelBottomUp(ch, getOrder(ch))));
}

}
