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
        public:
            Grid_1D() = default;
            Grid_1D(std::vector<node>& nodeV, std::vector<T1>& elements_,
                FiniteElement::Material::Material& mat_, std::vector<FiniteElement::boundary_node>& BCNodes,
                std::vector<FiniteElement::ExternalLoad>& ExternalLoads);
            ~Grid_1D() = default;

            auto Solve();

            auto GetLoadMatrix();
            auto GetBCMatrix();

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
            size_t m_maxBCNodes;
            Eigen::SparseMatrix<double> m_permColumn;
            Eigen::SparseMatrix<double> m_permRow;

            size_t dim;

            auto SetNodes(std::vector<node>& nodes_);
            auto SetElements(std::vector<T1>& elements_);
            auto SetMaterial(FiniteElement::Material::Material& material);
            auto SetBC(std::vector<FiniteElement::boundary_node>& BCNodes);
            auto SetLoads(std::vector<FiniteElement::ExternalLoad>& ExternalLoads);

            auto AssembleMatrices();
            auto AssembleLoadMatrix();
            auto AssemblePermutationMatrix();
            auto SolveSystem();
        };

        template<typename T1> Grid_1D<T1>::Grid_1D(std::vector<node> &nodeV, std::vector<T1>& elements_, 
            FiniteElement::Material::Material &mat_, std::vector<FiniteElement::boundary_node>& BCNodes, 
            std::vector<FiniteElement::ExternalLoad>& ExternalLoads)
        {
            dim = nodeV.size();
            m_material = new FiniteElement::Material::Material(mat_);
            m_MassMatrix.resize(dim, dim);
            m_MassMatrix.setZero();
            m_StiffMatrix.resize(dim, dim);
            m_StiffMatrix.setZero();
            m_LoadMatrix.resize(dim,1);
            m_LoadMatrix.setZero();
            m_permColumn.resize(dim, dim);
            m_permRow.resize(dim, dim);


            SetNodes(nodeV);
            SetElements(elements_);
            SetMaterial(mat_);
            SetBC(BCNodes);
            SetLoads(ExternalLoads);
        };

        template<typename T1> auto Grid_1D<T1>::SetNodes(std::vector<node>& nodes_)
        {
            nodeVector = &nodes_;
        };


        template<typename T1> auto Grid_1D<T1>::SetElements( std::vector<T1>& elements_)
        { 
            m_numElements = elements_.size();
            m_elements = &elements_;

            for (auto& i : (*m_elements))
            {
                i.SetNodeVector(nodeVector);
            }
        };

        template<typename T1> auto Grid_1D<T1>::SetMaterial(FiniteElement::Material::Material &material)
        { 
            for(size_t i=0; i < m_numElements; i++)
            {
                (*m_elements)[i].AttibuteMaterial(material);
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


        template<typename T1> auto Grid_1D<T1>::Solve()
        {

            AssembleMatrices();
            SolveSystem();

        }

        template<typename T1> auto Grid_1D<T1>::AssembleMatrices()
        {

            std::vector<size_t> IndexNodes;

            std::cout << "Mass Matrix - Before Assemble " << std::endl;
            std::cout << m_MassMatrix << std::endl;
            std::cout << "Stiff Matrix - Before Assemble " << std::endl;
            std::cout << m_StiffMatrix << std::endl;
            std::cout << "Load Matrix - Before Assemble " << std::endl;
            std::cout << m_LoadMatrix << std::endl;

            for (size_t i = 0; i < m_numElements; i++)
            {
                (*m_elements)[i].AssembleMassMatrix();
                (*m_elements)[i].AssembleStiffnessMatrix();

                auto tempMass = (*m_elements)[i].GetMassMatrix();
                auto tempStiff = (*m_elements)[i].GetStiffnessMatrix();
                IndexNodes = (*m_elements)[i].GetElementNodeIndex();

                auto numberNodesPerElement = IndexNodes.size();

                for (size_t j = 0; j < numberNodesPerElement; j++)
                {
                    auto jNode = IndexNodes[j];
                    for (size_t k = 0; k < numberNodesPerElement; k++)
                    {
                        auto kNode = IndexNodes[k];
                        //std::cout << tempMass(j * 2 + k) << std::endl;
                        m_MassMatrix(jNode, kNode) += tempMass(j * numberNodesPerElement + k);
                        m_StiffMatrix(jNode, kNode) += tempStiff(j * numberNodesPerElement + k);
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

            AssemblePermutationMatrix();

            m_MassMatrixReduced = m_permRow * m_MassMatrixReduced * m_permRow.transpose();
            m_StiffMatrixReduced = m_permRow * m_StiffMatrixReduced * m_permRow.transpose();

            std::cout << "Mass Matrix - Permuted" << std::endl;
            std::cout << m_MassMatrixReduced << std::endl;

            std::cout << "Stiff Matrix - Permuted" << std::endl;
            std::cout << m_StiffMatrixReduced << std::endl;


            std::cout << "Load Matrix - After Reduced " << std::endl;
            std::cout << m_LoadMatrixReduced << std::endl;

        };

        template<typename T1> auto Grid_1D<T1>::AssemblePermutationMatrix()
        {


            std::vector<size_t> indicesBCNodes;
            for (auto i : (*m_BcNodes))
            {
                if (i.GetBCType() == FiniteElement::bcType::dirichlet)
                {
                    indicesBCNodes.push_back(i.GetNode());
                };

            };

            std::vector<size_t> v((*m_BcNodes).size());
            std::vector<size_t>::iterator it;

            //sort indicesBCNodes
            std::sort(indicesBCNodes.begin(), indicesBCNodes.end());

            //it = std::set_union(indicesBCNodes.begin(), indicesBCNodes.end(), indicesBCNodes.begin(), indicesBCNodes.end(), v.begin());
            //v.resize(it - v.begin());

            m_maxBCNodes = indicesBCNodes.size();
            size_t controlDown = dim - m_maxBCNodes;
            size_t controlUp = 0;
            size_t control = 0;

            it = indicesBCNodes.begin();

            for (size_t i = 0; i < dim; i++) //iter in nodeVector
            {
                if (i != *it)
                {
                    m_permRow.coeffRef(controlUp, control) = 1.0;
                    m_permColumn.coeffRef(control, controlUp) = 1.0;
                    controlUp++;
                    control++;
                }
                else
                {
                    m_permRow.coeffRef(controlDown, i) = 1.0;
                    m_permColumn.coeffRef(i, controlDown) = 1.0;
                    controlDown++;
                    control++;

                    if (*it < (indicesBCNodes.size() - 1))
                    {
                        it++;
                    }
                }
            }

            std::cout << "Permutation Matrix" << std::endl;
            std::cout << m_permRow << std::endl;

            std::cout << "Permutation Matrix" << std::endl;
            std::cout << m_permColumn << std::endl;
        };

        template<typename T1> auto Grid_1D<T1>::SolveSystem()
        {

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(m_StiffMatrixReduced.block(0, 0, dim - m_maxBCNodes, dim - m_maxBCNodes));
            solver.factorize(m_StiffMatrixReduced.block(0, 0, dim - m_maxBCNodes, dim - m_maxBCNodes));
            auto x = solver.solve(m_LoadMatrixReduced.block(0, 0, dim - m_maxBCNodes, 1));

            std::cout << "Resposta??? " << std::endl;
            std::cout << x << std::endl;

        };

        template<typename T1> auto Grid_1D<T1>::AssembleLoadMatrix()
        {   

            std::vector<size_t> IndexNodes;
            std::vector<double> tempLoad;

            auto m_numLoadsElem = (*m_ExternalLoad).size();

            for(size_t i=0; i < m_numLoadsElem; i++)
            {
                auto elem = (*m_ExternalLoad)[i].GetElem();
                auto load = (*m_ExternalLoad)[i].GetExternalLoad();
                (*m_elements)[elem].AssembleLoadMatrix(load);
                auto tempLoad = (*m_elements)[elem].GetLoadMatrix();

                IndexNodes = (*m_elements)[elem].GetElementNodeIndex();

                for(size_t j=0; j < IndexNodes.size(); j++ )
                {
                    m_LoadMatrix(IndexNodes[j]) += tempLoad[j];
                }
                
            };

            auto m_numLoadsNode = (*m_BcNodes).size();

            for (size_t i = 0; i < m_numLoadsNode; i++)
            {
                if ((*m_BcNodes)[i].GetBCType() == FiniteElement::bcType::neumann)
                {
                    auto nodeLoad = (*m_BcNodes)[i].GetNode();
                    auto load = (*m_BcNodes)[i].GetBCValue();

                    m_LoadMatrix(nodeLoad) += load;
                }

            };

        };
        
    }
    
};