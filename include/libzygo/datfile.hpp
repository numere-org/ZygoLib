/*****************************************************************************
    ZygoLib
    Copyright (C) 2024  Erik Haenel et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef DATFILE_HPP
#define DATFILE_HPP

#include <string>
#include <vector>
#include <fstream>

#include "interferogramdata.hpp"

namespace ZygoLib
{
    /////////////////////////////////////////////////
    /// \brief Implementation of the ZYGO-Dat file
    /// type.
    /////////////////////////////////////////////////
    class DatFile
    {
        private:
            std::fstream datFileStream;

            uint16_t read2Bytes();
            uint32_t read4Bytes();
            DatFileHeader readHeader();
            IntMatrix readIntensityMatrix(const DatFileHeader& header);
            PhaseMatrix readPhaseMatrix(const DatFileHeader& header);

            double getScaleFactor(const DatFileHeader& header);

            void write2Bytes(uint16_t bytes);
            void write4Bytes(uint32_t bytes);
            void writeNBytes(size_t n);
            void writeHeader(const DatFileHeader& header);
            void writeIntensityMatrix(const IntMatrix& intMatrix);
            void writePhaseMatrix(const PhaseMatrix& phaseMatrix, double scale);

        public:
            DatFile(const std::string& sFileName, bool trunc = false);
            std::vector<InterferogramData> read();
            void write(const std::vector<InterferogramData>& vData);

            static bool isDatFile(const std::string& sFileName);
    };
}

#endif // DATFILE_HPP


