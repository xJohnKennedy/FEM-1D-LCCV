#include "BaseElement.hpp"


#pragma once


namespace  FiniteElement
{

    class OneDimensionalElement: public BaseElement
        {
        private:

        protected:

        public:
            OneDimensionalElement() = default;
            ~OneDimensionalElement() = default;
            OneDimensionalElement(std::vector<node> &nodeV, std::vector<size_t> IndexNodes);
            double GetLenght();
            Eigen::Array2d ShapeFunction(double ParametricPostition);
            Eigen::Array2d ShapeFunctionDerivative(double ParametricPostition);

        };
        
        OneDimensionalElement::OneDimensionalElement(std::vector<node> &nodeV, std::vector<size_t> IndexNodes) : BaseElement(nodeV,IndexNodes){};

        double OneDimensionalElement::GetLenght()
        {
            auto a0 = (*nodeVector)[index[0]];
            auto a1 = (*nodeVector)[index[index.size()-1]];
            return a1.GetPosition() - a0.GetPosition() ;
        };

        Eigen::Array2d OneDimensionalElement::ShapeFunction(double ParametricPostition = 0)
        {
            double l = GetLenght();
            return Eigen::Array2d { (1 - (ParametricPostition + 1) / 2) , ((ParametricPostition + 1) / 2) };
        };

        Eigen::Array2d OneDimensionalElement::ShapeFunctionDerivative(double ParametricPostition = 0)
        {
            double l = GetLenght();
            return Eigen::Array2d { -1.0f / l, 1.0f /  l};
        };
    
};