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

#include "include/libzygo/datfile.hpp"
#include <exception>
#include <cmath>

#define PHASENAN 0x7FFFFFF8

namespace ZygoLib
{
    /////////////////////////////////////////////////
    /// \brief Calculates the average of a vector
    /// that contains both NAN and numerical values.
    ///
    /// \param values const std::vector<double>&
    /// \return double
    ///
    /////////////////////////////////////////////////
    static double nanAvg(const std::vector<double>& values)
    {
        double sum = 0.0;
        size_t c = 0;

        for (double val : values)
        {
            if (!std::isnan(val))
            {
                sum += val;
                c++;
            }
        }

        if (c)
            return sum / c;

        return sum;
    }


    /////////////////////////////////////////////////
    /// \brief Advanced int cast
    ///
    /// \param number double
    /// \return int
    /// \remark Ported from NumeRe
    ///
    /////////////////////////////////////////////////
    static int intCast(double number)
    {
        // if quite close, use rint
        if (fabs(number - rint(number)) < 1e-7)
            return rint(number);

        // otherwise truncate
        return static_cast<int>(number);
    }


    /////////////////////////////////////////////////
    /// \brief Advanced matrix read (will return non-
    /// existent values as NaN)
    ///
    /// \param mat const std::vector<std::vector<T>>&
    /// \param row int
    /// \param col int
    /// \return double
    /// \remark Ported from NumeRe
    ///
    /////////////////////////////////////////////////
    template <class T>
    static double readMat(const std::vector<std::vector<T>>& mat, int row, int col)
    {
        if (row < (int)mat.size() && col < (int)mat[0].size() && row >= 0 && col >= 0)
            return mat[row][col];
        else
            return NAN;
    }


    /////////////////////////////////////////////////
    /// \brief Bilinear interpolation
    ///
    /// \param mat const std::vector<std::vector<T>>&
    /// \param row double
    /// \param col double
    /// \return double
    /// \remark Ported from NumeRe
    ///
    /////////////////////////////////////////////////
    template <class T>
    static double bilinearInterpol(const std::vector<std::vector<T>>& mat, double row, double col)
    {
        if (std::isnan(row) || std::isnan(col))
            return NAN;

        // Find the base index
        int nBaseLine = intCast(row) + (row < 0 ? -1 : 0);
        int nBaseCol = intCast(col) + (col < 0 ? -1 : 0);

        // Get the decimal part of the double indices
        double x = row - nBaseLine;
        double y = col - nBaseCol;

        // Find the surrounding four entries
        double f00 = readMat(mat, nBaseLine, nBaseCol);
        double f10 = readMat(mat, nBaseLine+1, nBaseCol);
        double f01 = readMat(mat, nBaseLine, nBaseCol+1);
        double f11 = readMat(mat, nBaseLine+1, nBaseCol+1);

        // If all are NAN, return NAN
        if (std::isnan(f00) && std::isnan(f01) && std::isnan(f10) && std::isnan(f11))
            return NAN;

        // Get average of non NAN values
        double avg = nanAvg({f00, f01, f10, f11});

        // Otherwise set NAN to average of all values
        f00 = std::isnan(f00) ? avg : f00;
        f10 = std::isnan(f10) ? avg : f10;
        f01 = std::isnan(f01) ? avg : f01;
        f11 = std::isnan(f11) ? avg : f11;

        //     f(0,0) (1-x) (1-y) + f(1,0) x (1-y) + f(0,1) (1-x) y + f(1,1) x y
        return f00*(1-x)*(1-y)    + f10*x*(1-y)    + f01*(1-x)*y    + f11*x*y;
    }


    /////////////////////////////////////////////////
    /// \brief Automatically rescale the moved matrix
    /// with the minimum of pixel height and width.
    ///
    /// \param matrix std::vector<std::vector<T>>&&
    /// \param pixel_height float
    /// \param pixel_width float
    /// \return std::vector<std::vector<T>>
    ///
    /////////////////////////////////////////////////
    template <class T>
    static std::vector<std::vector<T>> rescalePixels(std::vector<std::vector<T>>&& matrix, float pixel_height, float pixel_width)
    {
        if (pixel_height == pixel_width || !matrix.size())
            return matrix;

        float commonRes = std::min(pixel_width, pixel_height);
        float width_scale = pixel_width / commonRes;
        float height_scale = pixel_height / commonRes;
        size_t new_width = std::rint(matrix[0].size() * width_scale);
        size_t new_height = std::rint(matrix.size() * height_scale);

        std::vector<std::vector<T>> rescaled(new_height, std::vector<T>(new_width));

        // Resample the matrix
        for (size_t i = 0; i < rescaled.size(); i++)
        {
            for (size_t j = 0; j < rescaled[i].size(); j++)
            {
                // Find the (closest) source pixel by calculating the offset form the
                // center, scale it and add the center in the original image
                double source_row = i / height_scale;
                double source_col = j / width_scale;

                // Calculate the bilinear interpolation
                double newVal = bilinearInterpol(matrix, source_row, source_col);

                if (std::isnan(newVal) && !std::is_floating_point<T>::value)
                    rescaled[i][j] = 0;
                else
                    rescaled[i][j] = newVal;
            }
        }

        return rescaled;
    }


    /////////////////////////////////////////////////
    /// \brief DatFile constructor. Will open the
    /// target file on-the-fly.
    ///
    /// \param sFileName const std::string&
    /// \param trunc bool, true, when the output file
    /// shall be emptied first
    ///
    /////////////////////////////////////////////////
    DatFile::DatFile(const std::string& sFileName, bool trunc)
    {
        if (trunc)
            datFileStream.open(sFileName.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        else
            datFileStream.open(sFileName.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    }


    /////////////////////////////////////////////////
    /// \brief Reads two bytes from the stream and
    /// returns them as an unsigned 16-bit integer.
    ///
    /// \return uint16_t
    ///
    /////////////////////////////////////////////////
    uint16_t DatFile::read2Bytes()
    {
        uint16_t n2Bytes;
        datFileStream.read((char*)&n2Bytes, sizeof(n2Bytes));
        return __builtin_bswap16(n2Bytes);
    }


    /////////////////////////////////////////////////
    /// \brief Reads four bytes from the stream and
    /// returns them as an unsigned 32-bit integer.
    ///
    /// \return uint32_t
    ///
    /////////////////////////////////////////////////
    uint32_t DatFile::read4Bytes()
    {
        uint32_t n4Bytes;
        datFileStream.read((char*)&n4Bytes, sizeof(n4Bytes));
        return __builtin_bswap32(n4Bytes);
    }


    /////////////////////////////////////////////////
    /// \brief Reads the header of the Zygo-Dat-file.
    ///
    /// \return DatFileHeader
    ///
    /////////////////////////////////////////////////
    DatFileHeader DatFile::readHeader()
    {
        size_t curPos = datFileStream.tellg();

        uint32_t magic_number = read4Bytes();
        uint16_t header_format = read2Bytes();
        uint32_t header_size = read4Bytes();

        // Evaluate whether the header is correct
        // These are the only valid combinations
        // according to the file format documentation
        switch (header_format)
        {
            case 1:
                if (header_size != 834u || magic_number != 0x881B036F)
                    return DatFileHeader();
                break;
            case 2:
                if (header_size != 834u || magic_number != 0x881B0370)
                    return DatFileHeader();
                break;
            case 3:
                if (header_size != 4096u || magic_number != 0x881B0371)
                    return DatFileHeader();
                break;
            default:
                return DatFileHeader();
        }

        // Jump to first relevant header byte
        datFileStream.seekg(curPos + 48);

        DatFileHeader header;

        // Read intensity header part
        header.ac_org_x = read2Bytes();
        header.ac_org_y = read2Bytes();
        header.ac_width = read2Bytes();
        header.ac_height = read2Bytes();
        header.ac_n_buckets = read2Bytes();
        header.ac_range = read2Bytes();
        header.ac_n_bytes = read4Bytes();

        uint32_t n4Bytes;

        // Read phase header part
        header.cn_org_x = read2Bytes();
        header.cn_org_y = read2Bytes();
        header.cn_width = read2Bytes();
        header.cn_height = read2Bytes();
        header.cn_n_bytes = read4Bytes();
        datFileStream.seekg(curPos + 164);
        header.intf_scale_factor = *(float*)&(n4Bytes = read4Bytes());
        header.wavelength_in = *(float*)&(n4Bytes = read4Bytes());
        header.num_aperture = *(float*)&(n4Bytes = read4Bytes());
        header.obliquity_factor = *(float*)&(n4Bytes = read4Bytes());
        header.magnification = *(float*)&(n4Bytes = read4Bytes());
        header.camera_res = *(float*)&(n4Bytes = read4Bytes());
        datFileStream.seekg(curPos + 218);
        header.phase_res = read2Bytes();
        datFileStream.seekg(curPos + 360);
        header.refractive_index = *(float*)&(n4Bytes = read4Bytes());
        datFileStream.seekg(curPos + 486);
        header.pixel_width = *(float*)&(n4Bytes = read4Bytes());
        header.pixel_height = *(float*)&(n4Bytes = read4Bytes());
        datFileStream.seekg(curPos + 670);
        // LITTLE ENDIAN???? WHY???
        datFileStream.read((char*)&header.surface_type, sizeof(header.surface_type));

        // Jump to the first data matrix
        datFileStream.seekg(curPos + header_size);

        return header;
    }


    /////////////////////////////////////////////////
    /// \brief Reads the intensity matrix from the
    /// Zygo-Dat-file.
    ///
    /// \param header const DatFileHeader&
    /// \return IntMatrix
    /// \warning This matrix might be 3D
    ///
    /////////////////////////////////////////////////
    IntMatrix DatFile::readIntensityMatrix(const DatFileHeader& header)
    {
        IntMatrix intMatrix(header.ac_height, std::vector<uint16_t>(header.ac_width));

        for (uint16_t i = 0; i < header.ac_height; i++)
        {
            for (uint16_t j = 0; j < header.ac_width; j++)
            {
                intMatrix[i][j] = read2Bytes();
            }
        }

        return intMatrix;
    }


    /////////////////////////////////////////////////
    /// \brief Reads the phase matrix from the Zygo-
    /// Dat-file.
    ///
    /// \param header const DatFileHeader&
    /// \return PhaseMatrix
    ///
    /////////////////////////////////////////////////
    PhaseMatrix DatFile::readPhaseMatrix(const DatFileHeader& header)
    {
        PhaseMatrix phaseMatrix(header.cn_height, std::vector<double>(header.cn_width));

        double scale = getScaleFactor(header);

        for (uint16_t i = 0; i < header.cn_height; i++)
        {
            for (uint16_t j = 0; j < header.cn_width; j++)
            {
                uint32_t bytes = read4Bytes();

                if (bytes == PHASENAN)
                    phaseMatrix[i][j] = NAN;
                else
                    phaseMatrix[i][j] = (*(int32_t*)&bytes) * scale;
            }
        }

        return phaseMatrix;
    }


    /////////////////////////////////////////////////
    /// \brief Helper method to obtain the correct
    /// scaling factor.
    ///
    /// \param header const DatFileHeader&
    /// \return double
    ///
    /////////////////////////////////////////////////
    double DatFile::getScaleFactor(const DatFileHeader& header)
    {
        if (header.phase_res == 2)
            return header.intf_scale_factor * header.obliquity_factor * header.wavelength_in / 131072.0;
        else if (header.phase_res == 1)
            return header.intf_scale_factor * header.obliquity_factor * header.wavelength_in / 32768.0;
        else if (header.phase_res == 0) // Delphi code explicitly used 0
            return header.intf_scale_factor * header.obliquity_factor * header.wavelength_in / 4096.0;

        // In all other cases, the divisor seems to be 1 in Delphi
        return header.intf_scale_factor * header.obliquity_factor * header.wavelength_in;
    }


    /////////////////////////////////////////////////
    /// \brief Returns a Zygo-Dat-file and returns
    /// all available layers as a vector.
    ///
    /// \return std::vector<InterferogramData>
    ///
    /////////////////////////////////////////////////
    std::vector<InterferogramData> DatFile::read()
    {
        std::vector<InterferogramData> vData;

        while (!datFileStream.eof())
        {
            DatFileHeader rdHeader = readHeader();

            // Check for valid header values
            if (!rdHeader.ac_n_bytes && !rdHeader.cn_n_bytes)
                break;

            if (rdHeader.ac_n_buckets > 1)
                throw std::length_error("DatFile: Intensity matrix is 3D, which is not supported by ZygoLib.");

            // Create a new data set and read the
            // contents to memory
            vData.push_back(InterferogramData());

            vData.back().header = rdHeader;
            vData.back().intensityMatrix = readIntensityMatrix(rdHeader);
            vData.back().phaseMatrix = readPhaseMatrix(rdHeader);

            // Read and automatically transform the matrices to cope with different
            // pixel heights widths
            /*vData.back().intensityMatrix = rescalePixels(readIntensityMatrix(rdHeader), rdHeader.pixel_height, rdHeader.pixel_width);
            vData.back().phaseMatrix = rescalePixels(readPhaseMatrix(rdHeader), rdHeader.pixel_height, rdHeader.pixel_width);

            // Update the sizes correspondingly
            vData.back().header.setIntensitySizes(PixCoord(),
                                                  vData.back().intensityMatrix.size(), vData.back().intensityMatrix.front().size());
            vData.back().header.setPhaseSizes(PixCoord(),
                                              vData.back().phaseMatrix.size(), vData.back().phaseMatrix.front().size());

            // Ensure that the stored camera resolution actually reflects the
            // rescaled resolution
            if (rdHeader.pixel_height > 0 && rdHeader.pixel_width > 0)
                vData.back().header.camera_res = std::min(rdHeader.pixel_height, rdHeader.pixel_width);*/
        }

        return vData;
    }


    /////////////////////////////////////////////////
    /// \brief Writes the two bytes of the 16-bit
    /// integer to the output stream.
    ///
    /// \param bytes uint16_t
    /// \return void
    ///
    /////////////////////////////////////////////////
    void DatFile::write2Bytes(uint16_t bytes)
    {
        bytes = __builtin_bswap16(bytes);
        datFileStream.write((char*)&bytes, sizeof(bytes));
    }


    /////////////////////////////////////////////////
    /// \brief Writes the four bytes of the 32-bit
    /// integer to the output stream.
    ///
    /// \param bytes uint32_t
    /// \return void
    ///
    /////////////////////////////////////////////////
    void DatFile::write4Bytes(uint32_t bytes)
    {
        bytes = __builtin_bswap32(bytes);
        datFileStream.write((char*)&bytes, sizeof(bytes));
    }


    /////////////////////////////////////////////////
    /// \brief Writes a number of zero bytes. Those
    /// are needed to fill the header up to the
    /// desired size.
    ///
    /// \param n size_t
    /// \return void
    ///
    /////////////////////////////////////////////////
    void DatFile::writeNBytes(size_t n)
    {
        for (size_t i = 0; i < n; i++)
        {
            datFileStream.put(0u);
        }
    }


    /////////////////////////////////////////////////
    /// \brief Writes a header of the Zygo-Dat-file
    /// in version 1 format.
    ///
    /// \param header const DatFileHeader&
    /// \return void
    ///
    /////////////////////////////////////////////////
    void DatFile::writeHeader(const DatFileHeader& header)
    {
        // Write header signature
        write4Bytes(0x881B036F);
        write2Bytes(1u);
        write4Bytes(834u);

        // Fill up
        writeNBytes(38);

        // Write intensity part
        write2Bytes(header.ac_org_x);
        write2Bytes(header.ac_org_y);
        write2Bytes(header.ac_width);
        write2Bytes(header.ac_height);
        write2Bytes(header.ac_n_buckets);
        write2Bytes(header.ac_range);
        write4Bytes(header.ac_n_bytes);

        // Write phase part
        write2Bytes(header.cn_org_x);
        write2Bytes(header.cn_org_y);
        write2Bytes(header.cn_width);
        write2Bytes(header.cn_height);
        write4Bytes(header.cn_n_bytes);
        writeNBytes(164-76);
        write4Bytes(*(uint32_t*)&header.intf_scale_factor);
        write4Bytes(*(uint32_t*)&header.wavelength_in);
        write4Bytes(*(uint32_t*)&header.num_aperture);
        write4Bytes(*(uint32_t*)&header.obliquity_factor);
        write4Bytes(*(uint32_t*)&header.magnification);
        write4Bytes(*(uint32_t*)&header.camera_res);
        writeNBytes(218-188);
        write2Bytes(header.phase_res);
        writeNBytes(360-220);
        write4Bytes(*(uint32_t*)&header.refractive_index);
        writeNBytes(670-364);
        // LITTLE ENDIAN
        datFileStream.write((char*)&header.surface_type, sizeof(header.surface_type));
        writeNBytes(834-672);
    }


    /////////////////////////////////////////////////
    /// \brief Writes the intensity matrix to the
    /// output stream.
    ///
    /// \param intMatrix const IntMatrix&
    /// \return void
    /// \warning The intensity matrix might be 3D
    ///
    /////////////////////////////////////////////////
    void DatFile::writeIntensityMatrix(const IntMatrix& intMatrix)
    {
        for (size_t i = 0; i < intMatrix.size(); i++)
        {
            for (size_t j = 0; j < intMatrix[i].size(); j++)
            {
                write2Bytes(intMatrix[i][j]);
            }
        }
    }


    /////////////////////////////////////////////////
    /// \brief Writes the phase matrix to the output
    /// stream.
    ///
    /// \param phaseMatrix const PhaseMatrix&
    /// \param scale double
    /// \return void
    ///
    /////////////////////////////////////////////////
    void DatFile::writePhaseMatrix(const PhaseMatrix& phaseMatrix, double scale)
    {
        for (size_t i = 0; i < phaseMatrix.size(); i++)
        {
            for (size_t j = 0; j < phaseMatrix[i].size(); j++)
            {
                if (std::isnan(phaseMatrix[i][j]))
                    write4Bytes(PHASENAN);
                else
                {
                    int32_t bytes = (int32_t)(phaseMatrix[i][j] / scale);
                    write4Bytes(*(uint32_t*)&bytes);
                }
            }
        }
    }


    /////////////////////////////////////////////////
    /// \brief Writes a complete array of
    /// InterferogramData instances into a Zygo-Dat-
    /// file.
    ///
    /// \param vData const std::vector<InterferogramData>&
    /// \return void
    ///
    /////////////////////////////////////////////////
    void DatFile::write(const std::vector<InterferogramData>& vData)
    {
        for (size_t i = 0; i < vData.size(); i++)
        {
            writeHeader(vData[i].header);
            writeIntensityMatrix(vData[i].intensityMatrix);
            writePhaseMatrix(vData[i].phaseMatrix, getScaleFactor(vData[i].header));
        }
    }


    /////////////////////////////////////////////////
    /// \brief Determines, whether the passed file is
    /// a Zygo Datfile.
    ///
    /// \param sFileName const std::string&
    /// \return bool
    ///
    /////////////////////////////////////////////////
    bool DatFile::isDatFile(const std::string& sFileName)
    {
        DatFile file(sFileName, false);

        if (!file.datFileStream.good())
            return false;

        file.datFileStream.seekg(0, std::ios_base::end);
        size_t len = file.datFileStream.tellg();

        if (len)
        {
            file.datFileStream.seekg(0, std::ios_base::beg);

            uint32_t magic_number = file.read4Bytes();
            uint16_t header_format = file.read2Bytes();
            uint32_t header_size = file.read4Bytes();

            // Evaluate whether the header is correct
            // These are the only valid combinations
            // according to the file format documentation
            switch (header_format)
            {
                case 1:
                    if (header_size != 834u || magic_number != 0x881B036F)
                        return false;
                    break;

                case 2:
                    if (header_size != 834u || magic_number != 0x881B0370)
                        return false;
                    break;

                case 3:
                    if (header_size != 4096u || magic_number != 0x881B0371)
                        return false;
                    break;

                default:
                    return false;
            }

            return true;
        }

        return false;
    }
}


