#include "node.hpp"
#include "Element/element.hpp"
#include "load.hpp"
#include <memory>

#pragma once

namespace FiniteElement
{
    namespace Grid{

        template<typename T1> class Grid_1D
        {
        private:
            /* data */
            std::vector<node>* nodeVector;
            std::vector<FiniteElement::boundary_node>* m_BcNodes;
            std::vector<T1>* m_elements;
            int m_numElements;
            FiniteElement::Material::Material* m_material;
            Eigen::MatrixXd m_MassMatrix;
            Eigen::MatrixXd m_StiffMatrix;
            Eigen::MatrixXd m_LoadMatrix;
            Eigen::SparseMatrix<double> m_MassMatrixReduced;
            Eigen::SparseMatrix<double> m_StiffMatrixReduced;
            Eigen::MatrixXd m_LoadMatrixReduced;
            std::vector<FiniteElement::ExternalLoad>* m_ExternalLoad;
        public:
            Grid_1D() = default;
            Grid_1D(const std::vector<node> &nodeV, const FiniteElement::Material::Material &mat_);
            ~Grid_1D() = default;

            size_t dim;

            auto SetNodes(std::vector<node>& nodes_);
            auto SetElements( std::vector<T1>& elements_);
            auto SetMaterial(FiniteElement::Material::Material &material);
            auto SetBC(std::vector<FiniteElement::boundary_node> &BCNodes);
            auto SetLoads(std::vector<FiniteElement::ExternalLoad> &ExternalLoads);
            auto AssembleMesh();


            auto GetLoadMatrix();
            auto GetBCMatrix();
            auto AssembleLoadMatrix();

        };

        template<typename T1> Grid_1D<T1>::Grid_1D(const std::vector<node> &nodeV, const FiniteElement::Material::Material &mat_)
        {
            dim = nodeV.size();
            m_material = new FiniteElement::Material::Material(mat_);
            m_MassMatrix.resize(dim, dim);
            m_MassMatrix.setZero();
            m_StiffMatrix.resize(dim, dim);
            m_StiffMatrix.setZero();
            m_LoadMatrix.resize(dim,1);
            m_LoadMatrix.setZero();
        };

        template<typename T1> auto Grid_1D<T1>::SetNodes(std::vector<node>& nodes_)
        {
            nodeVector = &nodes_;
        };


        template<typename T1> auto Grid_1D<T1>::SetElements( std::vector<T1>& elements_)
        { 
            m_numElements = elements_.size();
            m_elements = &elements_;
        };

        template<typename T1> auto Grid_1D<T1>::SetMaterial(FiniteElement::Material::Material &material)
        { 
            for(size_t i=0; i < m_numElements; i++)
            {
                m_elements->operator[](i).AttibuteMaterial(material);
            };
        };

        
        template<typename T1> auto Grid_1D<T1>::SetBC(std::vector<FiniteElement::boundary_node> &BCNodes)
        { 
            m_BcNodes = &BCNodes;
        };

        template<typename T1> auto Grid_1D<T1>::SetLoads(std::vector<FiniteElement::ExternalLoad> &ExternalLoads)
        { 
            m_ExternalLoad = &ExternalLoads;
        };


        template<typename T1> auto Grid_1D<T1>::AssembleMesh()
        {   

            std::vector<size_t> IndexNodes;
            
            std::cout << "Mass Matrix - Before Assemble " << std::endl;
            std::cout << m_MassMatrix << std::endl;
            std::cout << "Stiff Matrix - Before Assemble " << std::endl;
            std::cout << m_StiffMatrix << std::endl;
            std::cout << "Load Matrix - Before Assemble " << std::endl;
            std::cout << m_LoadMatrix << std::endl;

            //TODO: create a method to assemble matrices
            for(size_t i=0; i < m_numElements; i++)
            {
                //TODO: after swapping pointers by references, adequate to call operator directly
                m_elements->operator[](i).AssembleMassMatrix();
                m_elements->operator[](i).AssembleStiffnessMatrix();

                auto tempMass = m_elements->operator[](i).GetMassMatrix();
                auto tempStiff = m_elements->operator[](i).GetStiffnessMatrix();
                IndexNodes = m_elements->operator[](i).GetElementNodeIndex();
                
                //TODO: rename sizeVector to numberNodesPerElement
                auto sizeVector = IndexNodes.size();

                for (size_t j = 0; j < sizeVector; j++)
                {
                    auto jNode = IndexNodes[j];
                    for (size_t k = 0; k < sizeVector; k++)
                    {
                        auto kNode = IndexNodes[k];
                        //std::cout << tempMass(j * 2 + k) << std::endl;
                        m_MassMatrix(jNode, kNode)     += tempMass(j* sizeVector + k);
                        m_StiffMatrix(jNode,kNode)    += tempStiff(j* sizeVector + k);
                    };
                };
            };

            std::cout << "Mass Matrix - After Assemble " << std::endl;
            std::cout << m_MassMatrix << std::endl;
            std::cout << "Stiff Matrix - After Assemble " << std::endl;
            std::cout << m_StiffMatrix << std::endl;


            this->AssembleLoadMatrix();

            std::cout << "Load Matrix - After Assemble " << std::endl;
            std::cout << m_LoadMatrix << std::endl;

            // convert dense matrix to sparse matrix
            m_MassMatrixReduced = m_MassMatrix.sparseView();
            m_StiffMatrixReduced = m_StiffMatrix.sparseView();
            m_LoadMatrixReduced = m_LoadMatrix;

            std::cout << "Mass Matrix - Before Reduced " << std::endl;
            std::cout << m_MassMatrixReduced << std::endl;



            //for (size_t temp = 0; temp < m_BcNodes->size(); temp++)
            //{
            //    auto nodeDirichlet = m_BcNodes->operator[](temp).GetNode();
            //    if (m_BcNodes->operator[](temp).GetBCType() == FiniteElement::bcType::dirichlet)
            //    {
            //        m_MassMatrixReduced.col(temp) *= 0;
            //        m_StiffMatrixReduced.col(temp) *= 0;
            //        m_StiffMatrixReduced.row(temp) *= 0;
            //        m_StiffMatrixReduced.coeffRef(1, 1);
            //        m_LoadMatrixReduced.row(temp) *= 0;
            //    }
            //};

            //code to permutation matrix
            //fist we need a indices of boundary nodes conditions

            //TODO: create a method to create permutation matrix
            Eigen::SparseMatrix<double> permColumn(dim, dim);
            Eigen::SparseMatrix<double> permRow(dim, dim);

            std::vector<size_t> indicesBCNodes;
            for (auto i : (*m_BcNodes))
            { 
                indicesBCNodes.push_back(i.GetNode());
            };
            
            //sort indicesBCNodes
            std::sort(indicesBCNodes.begin(), indicesBCNodes.end());
            size_t maxBCNodes = indicesBCNodes.size();
            size_t controlDown = dim - maxBCNodes;
            size_t controlUp = 0;
            size_t control = 0;

            std::vector<size_t>::iterator it = indicesBCNodes.begin();

            for (size_t i = 0; i < dim; i++) //iter in nodeVector
            {
                if (i != *it)
                {
                    permRow.coeffRef(controlUp, control) = 1.0;
                    permColumn.coeffRef(control, controlUp) = 1.0;
                    controlUp++;
                    control++;
                }
                else 
                {
                    permRow.coeffRef(controlDown, i) = 1.0;
                    permColumn.coeffRef(i, controlDown) = 1.0;
                    controlDown++;
                    control++;

                    if (*it < (indicesBCNodes.size() - 1))
                    {
                        it++;
                    }
                }


            }

            std::cout << "Permutation Matrix" << std::endl;
            std::cout << permRow << std::endl;

            std::cout << "Permutation Matrix" << std::endl;
            std::cout << permColumn << std::endl;


            m_MassMatrixReduced = permRow * m_MassMatrixReduced * permRow.transpose();
            m_StiffMatrixReduced = permRow * m_StiffMatrixReduced * permRow.transpose();

            std::cout << "Mass Matrix - Permuted" << std::endl;
            std::cout << m_MassMatrixReduced << std::endl;

            std::cout << "Stiff Matrix - Permuted" << std::endl;
            std::cout << m_StiffMatrixReduced << std::endl;


            std::cout << "Load Matrix - After Reduced " << std::endl;
            std::cout << m_LoadMatrixReduced << std::endl;

            //TODO: create a method to solve linear system
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(m_StiffMatrixReduced.block(0,0, dim - maxBCNodes, dim - maxBCNodes));
            solver.factorize(m_StiffMatrixReduced.block(0, 0, dim - maxBCNodes, dim - maxBCNodes));
            auto x = solver.solve(m_LoadMatrixReduced.block(0, 0, dim - maxBCNodes, 1));

            std::cout << "Resposta??? " << std::endl;
            std::cout << x << std::endl;

        };

        template<typename T1> auto Grid_1D<T1>::AssembleLoadMatrix()
        {   

            std::vector<size_t> IndexNodes;
            std::vector<double> tempLoad;

            auto m_numLoadsElem = m_ExternalLoad->size();

            for(size_t i=0; i < m_numLoadsElem; i++)
            {
                auto elem = m_ExternalLoad->operator[](i).GetElem();
                auto load = m_ExternalLoad->operator[](i).GetExternalLoad();
                m_elements->operator[](elem).AssembleLoadMatrix(load);
                auto tempLoad = m_elements->operator[](elem).GetLoadMatrix();

                IndexNodes = m_elements->operator[](elem).GetElementNodeIndex();

                for(size_t j=0; j < IndexNodes.size(); j++ )
                {
                    m_LoadMatrix(IndexNodes[j]) += tempLoad[j];
                }
                
            };

        };
        
    }
    
};