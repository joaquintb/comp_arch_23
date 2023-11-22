#include <gtest/gtest.h>

#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/simulation.hpp"
#include "../sim/progargs.hpp"
#include <cmath>
#include <fstream>
#include <sys/stat.h>


TEST(FileReadabilityTest, NonReadableFile) {
  // Create a temporary file with read-only permissions
  std::ofstream tmpFile("nonreadablefile.txt");
  tmpFile << "Non Readable";
  tmpFile.close();

  std::ifstream tmpInputFile("nonreadablefile.txt");
  tmpFile.close();

  //create temporary output file with write permissions
  std::ofstream tmpOutFile("tmpOutFile.txt");
  tmpOutFile << "Output File for Writing";
  tmpOutFile.close();


  // Make the file non-readable
  chmod("nonreadablefile.txt", 0222);  // Remove read permissions

  // Assert that the file is not readable, Need to fix inputfile to ifstream
  EXPECT_EXIT(check_files(tmpInputFile,"nonreadablefile.txt", tmpOutFile, "tmpOutFile.txt" ), testing::ExitedWithCode(-3), "Error: Cannot open nonreadablefile.txt for reading.\n");
}

//TEST(FileReadabilityTest, NonexistentFile) {
//  // Provide a file path that doesn't exist
//  std::string nonexistentFile = "nonexistentfile.txt";
//
//  // Assert that the nonexistent file is not readable
//  EXPECT_FALSE(isFileReadable(nonexistentFile));
//}
