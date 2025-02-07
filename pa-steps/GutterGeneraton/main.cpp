#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << argc << argv[0];
  std::cout << "GutterGeneration" << std::endl;
  std::ofstream out("GutterGeneration.txt");
  out << 1 << std::endl;
  out.close();
  return 0;
}