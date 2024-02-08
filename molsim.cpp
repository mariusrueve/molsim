#include <algorithm>
#include <iostream>
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include <exception>
#include <filesystem>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <spdlog/spdlog.h>

// Function to calculate the Tanimoto similarity between two molecules

int main(int argc, char *argv[]) {
    std::string inputPath;
    std::string databasePath;
    std::string outputPath;

    // create a description of the program options
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("input", boost::program_options::value<std::string>(&inputPath), "input file")
        ("database", boost::program_options::value<std::string>(&databasePath), "database file")
        ("output", boost::program_options::value<std::string>(&outputPath), "output file")
    ;

    // create a map to store the options
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm.count("input") == 0) {
        std::cerr << "Input file not specified" << std::endl;
        return 1;
    }

    if (vm.count("database") == 0) {
        std::cerr << "Database file not specified" << std::endl;
        return 1;
    }

    if (vm.count("output") == 0) {
        std::cerr << "Output file not specified" << std::endl;
        return 1;
    }

    // print out the options
    std::cout << "Input file: " << inputPath << std::endl;
    std::cout << "Database file: " << databasePath << std::endl;
    std::cout << "Output file: " << outputPath << std::endl;

    const auto inputFileExtension = std::filesystem::path(inputPath).extension().string();
    const auto databaseFileExtension = std::filesystem::path(databasePath).extension().string();

    RDKit::MOL_SPTR_VECT inputMols;
    RDKit::MOL_SPTR_VECT databaseMols;
    
    if (inputFileExtension == ".sdf") {
        RDKit::SDMolSupplier inputSupplier(inputPath, true, false, false);
        while (!inputSupplier.atEnd()) {
            auto mol = boost::make_shared<RDKit::RWMol>(*inputSupplier.next());
            if (mol) {
                inputMols.emplace_back(mol);
            }
        }
    } else {
        std::cerr << "Input file format not supported" << std::endl;
        return 1;
    }

    // print out the number of molecules in the input file
    std::cout << "Number of molecules in the input file: " << inputMols.size() << std::endl;

    return 0;
}