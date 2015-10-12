#include "evaluator.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>
#include <utility>

namespace io = boost::iostreams;
using namespace io;

namespace haplo {
    
double Evaluator::operator()(const char* ref_haplo_file, const char* sol_haplo_file) const 
{
    // Create some haplotypes
    haplo_container ref_haplotypes, sol_haplotypes;
    
    // Put some data into the haplotypes
    get_haplotypes(ref_haplo_file, ref_haplotypes);
    get_haplotypes(sol_haplo_file, sol_haplotypes);
    
    if (ref_haplotypes[0].size() != sol_haplotypes[0].size() &&
        ref_haplotypes[1].size() != sol_haplotypes[1].size() ) {
        std::cerr << "Different lengths!\n";
    }
    
    // Determine the score
    return 1.0 - recon_rate(ref_haplotypes, sol_haplotypes);
}

void Evaluator::get_haplotypes(const char* input_file, haplo_container& haplotypes) const 
{
    // Open file and convert to string (for tokenizer)
    io::mapped_file_source file(input_file);
    if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
    
    std::string data(file.data(), file.size());

    // Create a tokenizer to tokenize by newline character and by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> separator{"\n"};
    
    // Tokenize the data into lines
    tokenizer lines{data, separator};
    
    // Create a counter for the line number
    size_t counter = 0;
   
    // Get the data and store it in the data container 
    for (auto& line : lines) {
        // For each of the elements
        haplotypes[counter].resize(line.size());        
        for (auto i = 0; i < line.size(); ++i) {
            switch (line[i]) {
                case '0':
                    haplotypes[counter].set(i, 0);
                    break;
                case '1':
                    haplotypes[counter].set(i, 1);
                    break;
                default: 
                    std::cerr << "Error evaluating, incorrect input!\n";
                    break;
            }
        }
        ++counter;
    }
    
    if (file.is_open()) file.close();
}

double Evaluator::recon_rate(const haplo_container& ref_haplos, const haplo_container& sol_haplos) const 
{
    size_t count = 0;
    for (size_t i = 0; i < ref_haplos[0].size(); ++i) {
        count += std::min(  (ref_haplos[0].get(i) ^ sol_haplos[0].get(i)) + 
                            (ref_haplos[1].get(i) ^ sol_haplos[1].get(i))
                            ,
                            (ref_haplos[0].get(i) ^ sol_haplos[1].get(i)) + 
                            (ref_haplos[1].get(i) ^ sol_haplos[0].get(i)) 
                        );
    }
    return static_cast<double>(count) / ( 2.0 * static_cast<double>(ref_haplos[0].size()));
}

}       // End namespace haplo