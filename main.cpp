#include "src/FEM.hpp"

int main(int argc, char *argv[])
{
    // nodes
    //std::vector<FiniteElement::node> nodes{ 0.0, 1.0};
    std::vector<FiniteElement::node> nodes {0.0, 1.0, 2.0, 3.5};


    //elements type
    // std::vector <FiniteElement::Truss> barra1 {{nodes, {0,1}, 1.0f}};
    std::vector <FiniteElement::Truss> barra1{ {{0,1}, 2.0f}, {{1,2}, 1.0f}, {{1,2}, 1.0f}, {{1,3}, 1.0f}, {{2,3}, 1.0f} };

    //Material
    FiniteElement::Material::Homogeneous steel {7850, 210.0e10};

    //BC
    std::vector<FiniteElement::boundary_node> BCnodes { {0, FiniteElement::bcType::dirichlet, 0}, {0, FiniteElement::bcType::neumann, 10} };

    //External load
    std::vector<FiniteElement::ExternalLoad> Loads {{0, 10.f} };

    //mesh
    FiniteElement::Grid::Grid_1D<FiniteElement::Truss> mesh{nodes, barra1, steel, BCnodes, Loads};

    mesh.Solve();

    //linear solve
};