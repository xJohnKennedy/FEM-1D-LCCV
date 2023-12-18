#include "Element1D.hpp"
#include <iostream>



#pragma once

namespace  FiniteElement
{

    class Truss: public OneDimensionalElement 
    {
    private:
        double m_area = 0;
        Eigen::MatrixXd m_MassMatrix;
        Eigen::MatrixXd m_StiffMatrix;
        Eigen::MatrixXd m_LoadMatrix;
    public:
        Truss() = default;
        Truss(std::vector<size_t> IndexNodes, double area_);
        ~Truss() = default;


        double GetArea(double ParametricPostition);
        void AssembleMassMatrix();
        void AssembleStiffnessMatrix();
        void AssembleLoadMatrix (double &ExternalLoad);
        Eigen::ArrayXXd GetMassMatrix();
        Eigen::ArrayXXd GetStiffnessMatrix();
        Eigen::ArrayXd GetLoadMatrix();

    };    
    
};


FiniteElement::Truss::Truss(std::vector<size_t> IndexNodes, double area_): 
        FiniteElement::OneDimensionalElement(IndexNodes),
        m_area(area_){};

double FiniteElement::Truss::GetArea(double ParametricPostition = 0)
{
    return m_area;
}

void FiniteElement::Truss::AssembleMassMatrix ()
{
    Eigen::MatrixXd F0 = Truss::ShapeFunction(-0.57735);
    Eigen::MatrixXd F1 = Truss::ShapeFunction(+0.57735);
    double density0 = m_material->GetDensity(-0.57735);
    double density1 = m_material->GetDensity(+0.57735);
    double Area0 = GetArea(-0.57735);
    double Area1 = GetArea(+0.57735);

    std::cout << "Mass Matrix - Element " <<  std::endl;
    m_MassMatrix =   (F0 * F0.transpose() * density0 * Area0 +  F1 * F1.transpose() * density1 * Area1) * GetLenght() / 2 ;

    std::cout << m_MassMatrix << std::endl;
}

void FiniteElement::Truss::AssembleStiffnessMatrix ()
{
    Eigen::MatrixXd F0 = Truss::ShapeFunctionDerivative(-0.57735f);
    Eigen::MatrixXd F1 = Truss::ShapeFunctionDerivative(+0.57735f);
    std::vector<double> E0 = m_material->GetConstitutiveMatrix(-0.57735);
    std::vector<double> E1 = m_material->GetConstitutiveMatrix(+0.57735);
    double Area0 = GetArea(-0.57735);
    double Area1 = GetArea(+0.57735);
    
    std::cout << "Stiff Matrix - Element " << std::endl;
    m_StiffMatrix = (F0 * F0.transpose() * E0[0] * Area0 + F1 * F1.transpose() * E1[0] * Area1) * GetLenght() / 2;

    std::cout << m_StiffMatrix << std::endl;
}

void FiniteElement::Truss::AssembleLoadMatrix (double &ExternalLoad)
{
    Eigen::ArrayXd F0 = Truss::ShapeFunction(-0.57735f);
    Eigen::ArrayXd F1 = Truss::ShapeFunction(+0.57735f);

    std::cout << "Load Matrix - Element " << std::endl;
    m_LoadMatrix = Eigen::MatrixXd{   (F0 + F1)*ExternalLoad * GetLenght() / 2 };

    std::cout << m_LoadMatrix << std::endl;
}

Eigen::ArrayXXd FiniteElement::Truss::GetMassMatrix()
{
    return m_MassMatrix;
};


Eigen::ArrayXXd FiniteElement::Truss::GetStiffnessMatrix()
{
    return m_StiffMatrix;
};

Eigen::ArrayXd FiniteElement::Truss::GetLoadMatrix()
{
    return m_LoadMatrix;
};