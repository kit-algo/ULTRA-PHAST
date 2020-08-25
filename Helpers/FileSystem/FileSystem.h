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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <termios.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

#include <dirent.h>

#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/types.h>

#include "../Assert.h"
#include "../String/String.h"

namespace FileSystem {

    inline bool isAbsolutePath(const std::string& path) noexcept {
        return (!path.empty()) && (path[0] == '/');
    }

    inline bool isDirectory(const std::string& path) noexcept {
        DIR* dir;
        if ((dir = opendir(path.c_str())) != NULL) {
            closedir(dir);
            return true;
        } else {
            return false;
        }
    }

    inline bool isFile(const std::string& path) noexcept {
        std::ifstream f(path.c_str());
        if (f.good()) {
            f.close();
            return true;
        } else {
            return false;
        }
    }

    inline bool isFileOrDirectory(const std::string& path) noexcept {
        return isFile(path) || isDirectory(path);
    }

    inline std::string getParentDirectory(const std::string& fileName) noexcept {
        if (fileName.size() < 2) return "";
        size_t directoryEnd = fileName.find_last_of('/', fileName.size() - 2);
        if (directoryEnd >= fileName.size()) return "";
        return fileName.substr(0, directoryEnd);
    }

    inline std::string extendPath(const std::string& base, const std::string& file) noexcept {
        if (file.size() < 1) {
            return base;
        } else if (isAbsolutePath(file) || base.empty()) {
            return file;
        } else if (file[0] == '.') {
            if (file.size() < 2) {
                return base;
            } else if (file[1] == '/') {
                return extendPath(base, file.substr(2));
            }
            AssertMsg(file[1] == '.', "Cannot extend path '" << base << "' with file '" << file << "'!");
            if (file.size() == 2) {
                return getParentDirectory(base);
            }
            AssertMsg(file[2] == '/', "Cannot extend path '" << base << "' with file '" << file << "'!");
            return extendPath(getParentDirectory(base), file.substr(3));
        }
        if (base[base.length() - 1] == '/') {
            return base + file;
        } else {
            return base + "/" + file;
        }
    }

    inline std::string ensureExtension(const std::string& fileName, const std::string& extension) noexcept {
        if (String::endsWith(fileName, extension)) {
            return fileName;
        } else {
            return fileName + extension;
        }
    }

    inline std::string getFileNameWithoutExtension(const std::string& fileName) noexcept {
        size_t directoryEnd = fileName.find_last_of('/') + 1;
        if (directoryEnd >= fileName.size()) directoryEnd = 0;
        size_t extensionBegin  = fileName.find_last_of('.');
        if (extensionBegin <= directoryEnd) extensionBegin = fileName.size();
        return fileName.substr(directoryEnd, extensionBegin - directoryEnd);
    }

    inline void makeDirectory(const std::string& path) noexcept {
        if (path.empty()) return;
        if (isDirectory(path)) return;
        makeDirectory(getParentDirectory(path));
        mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    inline const std::string& ensureDirectoryExists(const std::string& fileName) noexcept {
        const std::string parentDirectory = getParentDirectory(fileName);
        if (parentDirectory == "") return fileName;
        if (isDirectory(parentDirectory)) return fileName;
        makeDirectory(parentDirectory);
        return fileName;
    }

}
