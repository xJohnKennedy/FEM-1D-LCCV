#include "src\FEM.hpp"

int main(int argc, char *argv[])
{
    // nodes
    std::vector<FiniteElement::node> nodes{ 0.0, 1.0};
    //std::vector<FiniteElement::node> nodes {0.0, 1.0, 2.0, 3.5};


    //elements type
    std::vector <FiniteElement::Truss> barra1 {{nodes, {0,1}, 1.0f}};
    //std::vector <FiniteElement::Truss> barra1{ {nodes, {0,1}, 2.0f}, {nodes, {1,2}, 1.0f}, {nodes, {1,2}, 1.0f}, {nodes, {1,3}, 1.0f}, {nodes, {2,3}, 1.0f} };

    //Material
    FiniteElement::Material::Homogeneous steel {7850, 210.0e10};
    
    //mesh
    FiniteElement::Grid::Grid_1D<FiniteElement::Truss> mesh{nodes, steel};

    //BC
    std::vector<FiniteElement::boundary_node> BCnodes { {0, FiniteElement::bcType::dirichlet, 0} };

    //External load
    std::vector<FiniteElement::ExternalLoad> Loads {{0, 10.f} };

    mesh.SetNodes(nodes);
    mesh.SetElements(barra1);
    mesh.SetMaterial(steel);
    mesh.SetBC(BCnodes);
    mesh.SetLoads(Loads);

    mesh.AssembleMesh();

    //linear solve


};