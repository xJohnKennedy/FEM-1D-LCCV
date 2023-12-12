#include "Element1D.hpp"


#pragma once

namespace  FiniteElement
{

    class Truss: public OneDimensionalElement 
    {
    private:
        double m_area = 0;
        std::vector<double> m_MassMatrix;
        std::vector<double> m_StiffMatrix;
    public:
        Truss() = default;
        Truss(std::vector<node> &nodeV, std::vector<size_t> IndexNodes, double area_);
        ~Truss() = default;


        double Area();
        void AssembleMassMatrix();
        void AssembleStiffnessMatrix();
        std::vector<double> GetMassMatrix();
        std::vector<double> GetStiffnessMatrix();

    };
    
};


FiniteElement::Truss::Truss(std::vector<node> &nodeV, std::vector<size_t> IndexNodes, double area_): 
        FiniteElement::OneDimensionalElement(nodeV, IndexNodes),
        m_area(area_){};

double FiniteElement::Truss::Area() 
{
    return m_area;
}

void FiniteElement::Truss::AssembleMassMatrix ()
{
    std::vector<double> F0 = Truss::ShapeFunction(-0.57735);
    std::vector<double> F1 = Truss::ShapeFunction(+0.57735);
    double density0 = m_material->GetDensity(-0.57735);
    double density1 = m_material->GetDensity(+0.57735);

    
    // copy elysion
    m_MassMatrix =  std::vector<double> {   F0[0]*F0[0]*density0 + F1[0]*F1[0]*density1,
                                            F0[0]*F0[1]*density0 + F1[0]*F1[1]*density1,
                                            F0[0]*F0[1]*density0 + F1[0]*F1[1]*density1,
                                            F0[1]*F0[1]*density0 + F1[1]*F1[1]*density1
                                        };
}

void FiniteElement::Truss::AssembleStiffnessMatrix ()
{
    std::vector<double> F0 = Truss::ShapeFunctionDerivative(-0.57735f);
    std::vector<double> F1 = Truss::ShapeFunctionDerivative(+0.57735f);
    std::vector<double> E0 = m_material->GetConstitutiveMatrix(-0.57735);
    std::vector<double> E1 = m_material->GetConstitutiveMatrix(+0.57735);


    m_StiffMatrix = std::vector<double> {   F0[0]*F0[0]*E0[0] + F1[0]*F1[0]*E1[0],
                                            F0[0]*F0[1]*E0[0] + F1[0]*F1[1]*E1[0],
                                            F0[0]*F0[1]*E0[0] + F1[0]*F1[1]*E1[0],
                                            F0[1]*F0[1]*E0[0] + F1[1]*F1[1]*E1[0]
                                        };
}

std::vector<double> FiniteElement::Truss::GetMassMatrix()
{
    return m_MassMatrix;
};


std::vector<double> FiniteElement::Truss::GetStiffnessMatrix()
{
    return m_StiffMatrix;
};