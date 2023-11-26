#pragma once

#include <vector>
#include <memory>
#include "../material.hpp"
#include "../node.hpp"

#include <Eigen\Sparse>

namespace  FiniteElement
{

    class BaseElement
    {

        // somentes nós
    private:
        
    protected:
        std::vector<node>* nodeVector;
        std::vector<size_t> index;
        FiniteElement::Material::Material* m_material;

    public:
        BaseElement() = default;
        ~BaseElement() = default;
        BaseElement(std::vector<node> &nodeV, std::vector<size_t> IndexNodes);

        std::vector<size_t> GetElementNodeIndex();
        
        virtual Eigen::Array2d ShapeFunction(double ParametricPostition) = 0;              // shape function
        virtual Eigen::Array2d ShapeFunctionDerivative(double ParametricPostition) = 0;    // shape function
        virtual  void AttibuteMaterial( FiniteElement::Material::Material &mat_);

    };
}

FiniteElement::BaseElement::BaseElement(std::vector<node> &nodeV, std::vector<size_t> IndexNodes = {}) : nodeVector(&nodeV), index(IndexNodes) {};

void FiniteElement::BaseElement::AttibuteMaterial( FiniteElement::Material::Material &mat_)
{
    m_material = &mat_;
};

std::vector<size_t> FiniteElement::BaseElement::GetElementNodeIndex()
{
    return index; 
};
