#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << argc << argv[0];
  std::cout << "Dissection" << std::endl;
  std::ofstream out("Dissection.txt");
  out << 1 << std::endl;
  out.close();
  return 0;
}