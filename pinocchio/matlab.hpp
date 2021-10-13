#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>
#include <fstream>

namespace
{
    const std::string name_tag = "# name: " ;
    const std::string type_tag = "# type: matrix" ;
    const std::string rows_tag = "# rows: 1" ;
}

// Save vector as matrix type in Octave (.mat) stm.
template <typename T>
std::ostream& save_vector_as_matrix( const std::string& name, const std::vector<T>& matrix, std::ofstream& stm )
{
    std::copy( matrix.begin(), matrix.end(), std::ostream_iterator<T>( stm, ", " ) ) ;
	return stm << "\n\n\n" ;
}

template <typename T>
std::istream& load_matrix_to_vector( std::vector<T>& matrix, const std::string& name, std::istream& stm )
{
    std::string str ;
    if ( !( std::getline( stm, str ) && str == ( name_tag + name ) ) ) goto failure ;
    if( !( std::getline( stm, str ) && str == type_tag ) ) goto failure ;
    if ( !( std::getline( stm, str ) && str == rows_tag  ) ) goto failure ;

    char ch;
    std::size_t expected_size ;
    if( !( stm >> ch && ch == '#' && stm >> str && str == "columns:" && stm >> expected_size ) )
        goto failure ;

    matrix.clear() ;
    using iterator = std::istream_iterator<T> ;
    std::copy( iterator(stm), iterator(), std::back_inserter(matrix) ) ;
    if( matrix.size() != expected_size ) goto failure ;

    return stm ;

    failure:
        stm.setstate( std::ios::failbit ) ;
        matrix.clear() ;
        return stm ;
}

template <typename T>
std::vector<T> load_matrix( const std::string& name, std::istream& stm )
{
    std::vector<T> matrix ;
    if( !load_matrix_to_vector( matrix, name, stm ).eof() ) std::cerr << "input failure!\n" ;
    return matrix ;
}