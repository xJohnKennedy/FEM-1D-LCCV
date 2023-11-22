#include "node.hpp"
#include "Element\element.hpp"
#include <memory>

#pragma once

namespace FiniteElement
{
    namespace Grid{

        template<typename T1> class Grid_1D
        {
        private:
            /* data */
            FiniteElement::boundary_node* BC_Nodes;
            T1* m_elements;
            int m_numElements;
            FiniteElement::Material::Material* m_material;
            double* m_MassMatrix;
            double* m_StiffMatrix;
        public:
            Grid_1D() = default;
            Grid_1D(const int &numElements_, const std::vector<node> &nodeV, const FiniteElement::Material::Material &mat_);
            ~Grid_1D() = default;

            size_t dim;

            auto SetElements(const T1* elements_);
            auto SetMaterial(FiniteElement::Material::Material* material);
            auto AssembleMesh();

        };

        template<typename T1> Grid_1D<T1>::Grid_1D(const int &numElements_, const std::vector<node> &nodeV, const FiniteElement::Material::Material &mat_)
        {
            dim = nodeV.size();
            m_numElements = numElements_;
            m_elements = new T1 [m_numElements];
            BC_Nodes = new FiniteElement::boundary_node[m_numElements + 1];
            m_material = new FiniteElement::Material::Material(mat_);
            m_MassMatrix = new double[dim*dim]();
            m_StiffMatrix = new double[dim*dim]();
        };


        template<typename T1> auto Grid_1D<T1>::SetElements(const T1* elements_)
        { 
            for(size_t i=0; i < m_numElements; i++)
            {
                m_elements[i] = elements_[i];
            };
        };

        template<typename T1> auto Grid_1D<T1>::SetMaterial(FiniteElement::Material::Material* material)
        { 
            for(size_t i=0; i < m_numElements; i++)
            {
                m_elements[i].AttibuteMaterial(material);
            };
        };

        template<typename T1> auto Grid_1D<T1>::AssembleMesh()
        {   

            std::vector<size_t> IndexNodes;
            std::vector<double> tempMass;
            std::vector<double> tempStiff;

            for(size_t i=0; i < m_numElements; i++)
            {
                m_elements[i].AssembleMassMatrix();
                m_elements[i].AssembleStiffnessMatrix();

                tempMass = m_elements[i].GetMassMatrix();
                tempStiff = m_elements[i].GetStiffnessMatrix();
                IndexNodes = m_elements[i].GetElementNodeIndex();

                for (size_t j = 0; j < 2; j++)
                {
                    auto jNode = IndexNodes[j];
                    for (size_t k = 0; k < 2; k++)
                    {
                        auto kNode = IndexNodes[k];
                        m_MassMatrix[jNode*dim + kNode]     += tempMass[j*2 + k];
                        m_StiffMatrix[jNode*dim + kNode]    += tempStiff[j*2 + k];
                    };
                };
            };
        };
        
    }
    
};