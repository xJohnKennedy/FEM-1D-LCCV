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
            std::vector<double> ShapeFunction(double ParametricPostition);
            std::vector<double> ShapeFunctionDerivative(double ParametricPostition);

        };
        
        OneDimensionalElement::OneDimensionalElement(std::vector<node> &nodeV, std::vector<size_t> IndexNodes) : BaseElement(nodeV,IndexNodes){};

        double OneDimensionalElement::GetLenght()
        {
            auto a0 = (*nodeVector)[index[0]];
            auto a1 = (*nodeVector)[index[index.size()-1]];
            return a1.GetPosition() - a0.GetPosition() ;
        };

        std::vector<double> OneDimensionalElement::ShapeFunction(double ParametricPostition = 0)
        {
            double l = GetLenght();
            return std::vector<double> { 1 - (ParametricPostition + 1) / 2 / l, (ParametricPostition + 1) / 2 / l};
        };

        std::vector<double> OneDimensionalElement::ShapeFunctionDerivative(double ParametricPostition = 0)
        {
            double l = GetLenght();
            return std::vector<double> { -1.0f / l, 1.0f /  l};
        };
    
};