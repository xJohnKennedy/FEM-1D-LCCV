#include "src\FEM.hpp"

int main(int argc, char *argv[])
{
    // nodes
    std::vector<FiniteElement::node> nodes {0.0, 1.0, 2.0};
    
    std::vector<size_t> index {0,1};

    //elements type
    FiniteElement::Truss barra1[2] {{nodes, {0,1}, 1.0f}, {nodes, {1,2}, 1.0f}};

    //Material
    FiniteElement::Material::Homogeneous steel {7850, 210.0e10};
    
    //mesh
    FiniteElement::Grid::Grid_1D<FiniteElement::Truss> mesh{2, nodes, steel};

    mesh.SetElements(&barra1[0]);
    mesh.SetMaterial(&steel);

    mesh.AssembleMesh();

};