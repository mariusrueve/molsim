#include <algorithm>
#include <iostream>
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include <exception>
#include <filesystem>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/SparseBitVect.h>
#include <DataStructs/BitVectUtils.h>

#include <spdlog/spdlog.h>

// Function to calculate the Tanimoto similarity between two molecules

int main(int argc, char *argv[])
{
    std::string inputPath;
    std::string databasePath;
    std::string outputPath;

    // create a description of the program options
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("input", boost::program_options::value<std::string>(&inputPath), "input file. Supported formats: .sdf")("database", boost::program_options::value<std::string>(&databasePath), "database file. Supported formats: .sdf");

    // create a map to store the options
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm.count("input") == 0)
    {
        std::cerr << "Input file not specified" << std::endl;
        return 1;
    }

    if (vm.count("database") == 0)
    {
        std::cerr << "Database file not specified" << std::endl;
        return 1;
    }

    // print out the options
    spdlog::info("Input file: {}", inputPath);
    spdlog::info("Database file: {}", databasePath);

    const auto inputFilename = std::filesystem::path(inputPath).filename().string();
    const auto databaseFilename = std::filesystem::path(databasePath).filename().string();

    const auto inputFileExtension = std::filesystem::path(inputPath).extension().string();
    const auto databaseFileExtension = std::filesystem::path(databasePath).extension().string();

    RDKit::MOL_SPTR_VECT inputMols;
    RDKit::MOL_SPTR_VECT databaseMols;

    if (inputFileExtension == ".sdf")
    {
        RDKit::SDMolSupplier inputSupplier(inputPath, true, false, false);
        while (!inputSupplier.atEnd())
        {
            auto mol = boost::make_shared<RDKit::RWMol>(*inputSupplier.next());
            if (mol)
            {
                inputMols.emplace_back(mol);
            }
        }
    }
    else if (inputFileExtension == ".smi")
    {
        RDKit::SmilesMolSupplier inputSupplier(inputPath, "\t", 0, 1, 0);
        while (!inputSupplier.atEnd())
        {
            auto mol = boost::make_shared<RDKit::RWMol>(*inputSupplier.next());
            if (mol)
            {
                inputMols.emplace_back(mol);
            }
        }
    }
    else
    {
        spdlog::error("Input file format not supported");
        return 1;
    }

    if (databaseFileExtension == ".sdf")
    {
        RDKit::SDMolSupplier databaseSupplier(databasePath, true, false, false);
        while (!databaseSupplier.atEnd())
        {
            auto mol = boost::make_shared<RDKit::RWMol>(*databaseSupplier.next());
            if (mol)
            {
                databaseMols.emplace_back(mol);
            }
        }
    }
    else if (databaseFileExtension == ".smi")
    {
        RDKit::SmilesMolSupplier databaseSupplier(databasePath, "\t", 0, 1, 0);
        while (!databaseSupplier.atEnd())
        {
            auto mol = boost::make_shared<RDKit::RWMol>(*databaseSupplier.next());
            if (mol)
            {
                databaseMols.emplace_back(mol);
            }
        }
    }
    else
    {
        spdlog::error("Database file format not supported");
        return 1;
    }

    // print out the number of molecules in the input file
    spdlog::info("Number of molecules in the input file: {}", inputMols.size());
    spdlog::info("Number of molecules in the database file: {}", databaseMols.size());

    std::unique_ptr<RDKit::FingerprintGenerator<std::uint32_t>> fpgen{RDKit::MorganFingerprint::getMorganGenerator<std::uint32_t>(2)};

    std::vector<ExplicitBitVect *> inputFingerprints;
    std::vector<ExplicitBitVect *> databaseFingerprints;

    for (const auto mol : inputMols)
    {
        inputFingerprints.emplace_back(fpgen->getFingerprint(*mol));
    }

    for (const auto mol : databaseMols)
    {
        databaseFingerprints.emplace_back(fpgen->getFingerprint(*mol));
    }

    spdlog::info("Finished generating all fingerprints.");

    std::vector<std::vector<double>> similarityMatrix(inputMols.size(), std::vector<double>(databaseMols.size(), 0.0));

    for (std::size_t i = 0; i < inputFingerprints.size(); ++i)
    {
        for (std::size_t j = 0; j < databaseFingerprints.size(); ++j)
        {
            similarityMatrix[i][j] = TanimotoSimilarity(*inputFingerprints[i], *databaseFingerprints[j]);
        }
    }

    spdlog::info("Finished calculating similarity matrix");

    // get the best match with the highest similarity for each input molecule
    std::vector<std::pair<std::size_t, std::size_t>> bestMatches;
    for (std::size_t i = 0; i < inputMols.size(); ++i)
    {
        auto maxSimilarity = std::max_element(similarityMatrix[i].begin(), similarityMatrix[i].end());
        auto index = std::distance(similarityMatrix[i].begin(), maxSimilarity);
        bestMatches.emplace_back(std::make_pair(i, index));
    }

    spdlog::info("Finished getting best matches");

    // write the best matches to a csv file in the format "smiles_input, smiles_database, similarity" with header "inputFilename, databaseFilename, similarity"
    std::ofstream outputFile("best_matches.csv");
    outputFile << inputFilename << ", " << databaseFilename << ", similarity" << std::endl;    
    for (const auto &match : bestMatches)
    {
        std::string smilesInput = RDKit::MolToSmiles(*inputMols[match.first]);
        std::string smilesDatabase = RDKit::MolToSmiles(*databaseMols[match.second]);
        outputFile << smilesInput << ", " << smilesDatabase << ", " << similarityMatrix[match.first][match.second] << std::endl;
    }
    outputFile.close();

    spdlog::info("Finished writing best matches to file: best_matches.csv");

    return 0;
}