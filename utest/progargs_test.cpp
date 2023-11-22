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
  // std::ofstream tmpFile("nonreadablefile.txt");
  // tmpFile << "Non Readable";
  // tmpFile.close();

  // std::ifstream tmpInputFile("nonreadablefile.txt");
  // tmpFile.close();

  const std::string filename = "nonexistent_file.txt";
  std::ifstream tmpInputFile(filename);
  tmpInputFile.close();


  //create temporary output file with write permissions
  std::ofstream tmpOutFile("tmpOutFile.txt");
  tmpOutFile << "Output File for Writing";
  tmpOutFile.close();


  // Make the file non-readable
  //chmod("nonreadablefile.txt", 222);  // Remove read permissions

  // Assert that the file is not readable, Need to fix inputfile to ifstream
  //need to use corresponding unsigned value, ex: -3 signed integer is 253 signed integer
  EXPECT_EXIT(check_files(tmpInputFile,"nonexistent_file.txt", tmpOutFile, "tmpOutFile.txt" ), testing::ExitedWithCode(253), "Error: Cannot open nonexistent_file.txt for reading.\n");
}

//TEST(FileReadabilityTest, NonexistentFile) {
//  // Provide a file path that doesn't exist
//  std::string nonexistentFile = "nonexistentfile.txt";
//
//  // Assert that the nonexistent file is not readable
//  EXPECT_FALSE(isFileReadable(nonexistentFile));
//}
