#include "BiljectivePara.h"
#include <fstream>
#include <io.h>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << "Syntax: " << argv[0] << " <input mesh>" << endl;
    return -1;
  }
  if (2 == argc) {
    const string input_mesh = argv[1];

    Mesh mesh;
    OpenMesh::IO::read_mesh(mesh, input_mesh);
    BiljectivePara *bil_para = new BiljectivePara(mesh, input_mesh);

    bil_para->load(); // 读入网格、完成shell的构建
    bil_para->parameterization();
    bil_para->write_obj("uv_mesh.obj");

    delete bil_para;
    bil_para = NULL;

    // system("pause");
    return 0;
  }
}
