/*  H5Tools
    defines standard procedures to save data in hdf5-files
    v0.1, 13.5.2020

 

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

*/
#include "h5tools.h"

H5::DataSet h5CreateDSet(H5::H5File *file, std::vector<hsize_t> dim,
                   std::string name, std::string type){
    /* create a chunked dataset (in order to enable compression)
     * the chunk size is the entire dataset 
     * Notes for compression:
     * Gzip:     deflation level in [0,9],  smaller is faster but less compressed
     * Szip:     1st. Argument = method selection 
     *              H5_SZIP_EC_OPTION_MASK = Entropy coding method
     *              OR
     *              H5_SZIP_NN_OPTION_MASK = Nearest neighbor mehtod
     *           2nd. Argument = pixel per block
     *               MUST be EVEN
     *               <=32
     *               <= # of elements in chunk
     *               best performance recommondation multiples of fastest changing dimension
     * Performance check: check compression of file with h5dump -pH ....h5
     */

    // enable chunking (dataset can be at minimum extended by chunk)
    // chunk size may not be too large to not cause an cache overflow
    H5::DSetCreatPropList cparms;
    std::vector<hsize_t> chunks;
    chunks = dim;
    int elems = 1;
    int divider = 2;
    for (unsigned int i=0; i<chunks.size(); i++)
        elems *= chunks[i];
    // std::cout << "chunk-elems: " << elems << "\n";
    while (elems > 5000){
        chunks[0] /= divider;       // reduce chunk in time-dimension
        elems /= divider;
        if (chunks[0] == 1)
            break;
    }
    // std::cout << "chunk-time: " << chunks[0] << "\n";
    
    cparms.setChunk(dim.size(), &chunks[0]);

    // select compression:
    cparms.setDeflate(9);        // Gzip compression
    // int PixelPerBlock = dim[dim.size()-1];
    // if (PixelPerBlock > 32)
    //         PixelPerBlock = 32;
    // if (PixelPerBlock < 32)
    //         PixelPerBlock *= (int)(32/PixelPerBlock);
    // if (PixelPerBlock%2 != 0)   // must be even
    //         PixelPerBlock += 1;
    // cparms.setSzip(H5_SZIP_EC_OPTION_MASK, PixelPerBlock);   
    // std::cout << "PixelPerBlock: " << PixelPerBlock << "\n";

	// Create the data space for the dataset.
    H5::DataSpace dataspace(dim.size(), &dim[0]);
    // create dataset 
    if (type == "int")
        return file->createDataSet(name.c_str(), H5::PredType::NATIVE_INT, dataspace, cparms);
    else        // if (type == "double")
        return file->createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace, cparms);
}

H5::DataSet h5CreateDSetExt(H5::H5File *file, std::vector<hsize_t> dim,
                     std::string name, std::string type){
    /* create and dataset which is extendable in the first dimension and 
     * whose chunk dimensions are the other dimensions left 
     * Notes for compression:
     * Gzip:     deflation level in [0,9],  smaller is faster but less compressed
     * Szip:     1st. Argument = method selection 
     *              H5_SZIP_EC_OPTION_MASK = Entropy coding method
     *              OR
     *              H5_SZIP_NN_OPTION_MASK = Nearest neighbor mehtod
     *           2nd. Argument = pixel per block
     *               <=32
     *               <= # of elements in chunk
     *               best performance recommondation multiples of fastest changing dimension
     * Performance check: check compression of file with h5dump -pH ....h5
     */

	// Create the data space for the dataset.
    std::vector<hsize_t> maxdim;
    maxdim = dim; maxdim[0] = H5S_UNLIMITED;
    H5::DataSpace dataspace(dim.size(), &dim[0], &maxdim[0]);

    // enable chunking (dataset can be at minimum extended by chunk)
    H5::DSetCreatPropList cparms;
    dim[0] = 1;
    cparms.setChunk(dim.size(), &dim[0]);

    // select compression:
    cparms.setDeflate(9);        // Gzip compression
    // cparms.setSzip(H5_SZIP_NN_OPTION_MASK, 20);   

    // create dataset 
    if (type == "int")
        return file->createDataSet(name.c_str(), H5::PredType::NATIVE_INT, dataspace, cparms);
    else        // if (type == "double")
        return file->createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace, cparms);
}

void h5WriteInt(H5::DataSet *dset, std::vector<int> data, std::vector<hsize_t> offset){

    // define dimension of data
    // ASSUME: always written to last dimension
    std::vector<hsize_t> dimdata(offset.size(), 1);
    dimdata[dimdata.size()-1] = data.size();

    // Select a hyperslab to copy chunk to:
    H5::DataSpace filespace = dset->getSpace ();
    filespace.selectHyperslab(H5S_SELECT_SET, &dimdata[0], &offset[0]); 

    // define memory of data
    H5::DataSpace mspace(offset.size(), &dimdata[0]);

    // write the data to dataset
    // 1st. filespace defines memoryspace of newdata and 2nd the filespace of the dataset
    dset->write(&data[0], H5::PredType::NATIVE_INT, mspace, filespace);

    filespace.close();
    mspace.close();
}

void h5WriteDouble(H5::DataSet *dset, std::vector<double> data, std::vector<hsize_t> offset){

    // define dimension of data
    // ASSUME: always written to last dimension
    std::vector<hsize_t> dimdata(offset.size(), 1);
    dimdata[dimdata.size()-1] = data.size();

    // Select a hyperslab to copy chunk to:
    H5::DataSpace filespace = dset->getSpace ();
    filespace.selectHyperslab(H5S_SELECT_SET, &dimdata[0], &offset[0]); 

    // define memory of data
    H5::DataSpace mspace(offset.size(), &dimdata[0]);

    // write the data to dataset
    // 1st. filespace defines memoryspace of newdata and 2nd the filespace of the dataset
    dset->write(&data[0], H5::PredType::NATIVE_DOUBLE, mspace, filespace);

    filespace.close();
    mspace.close();
}


void h5readDimension(H5::DataSet *dataset, std::vector<hsize_t> &dim){
    dim.resize(0);
    H5::DataSpace filespace = dataset->getSpace ();
    int ndims = filespace.getSimpleExtentNdims();

    dim.resize(ndims);
    filespace.getSimpleExtentDims(&dim[0]);
}

void h5read_extend_dim(H5::DataSet *dataset, std::vector<hsize_t> &dim){
    h5readDimension(dataset, dim);
    dim[0] += 1;   // dummy-prey
    dataset->extend(&dim[0]);
}

// takes path to f_h5file as input
void h5CreateWriteDset(std::string f_h5out, std::string n_group,
                       std::string n_dset,
                       std::vector< std::vector<double> > data,
                       bool extend){
    H5::H5File *H5out;
    H5::DataSet *dset;
    std::string n_full_dset;
    if (n_group == "xx")
        n_full_dset = "/" + n_dset;
    else
        n_full_dset = "/" + n_group + "/" + n_dset;
    // if ( boost::filesystem::exists(f_h5out) )
    if ( exists(f_h5out) )
        H5out = new H5::H5File(f_h5out.c_str(), H5F_ACC_RDWR);
    else
        H5out = new H5::H5File(f_h5out.c_str(), H5F_ACC_TRUNC);
    std::vector<hsize_t> dim(2, 0);
    if ( extend ){    // multiple Runs -> open and extend existing files
        dset = new H5::DataSet(H5out->openDataSet(n_full_dset.c_str()));
        h5read_extend_dim(dset, dim);
    }
    else{   // single run save the dataset as it is
        dim[0] = data.size();
        dim[1] = data[0].size();
        dset = new H5::DataSet(h5CreateDSet(H5out, dim,
                                            n_full_dset.c_str(), "double"));
    }
    std::vector<hsize_t> offset(dim.size());  // offset of hyperslab  
    if ( extend )
        offset[0] = dim[0]-1;  // to not overwrite preceding run
    for(int i=0; i<dim[dim.size()-2]; i++){
        h5WriteDouble(dset, data[i], offset);
        offset[offset.size()-2]++;
    }
    delete dset;
    delete H5out;
}
