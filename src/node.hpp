#pragma once

namespace FiniteElement {
        
    class node
    {
    private:
        double m_xPostion = 0 ;
    public:
        node() = default;
        node(const double &x_);
        ~node() = default;
        auto GetPosition();
    };

    node::node (const double &x_) : m_xPostion(x_) {};

    auto node::GetPosition()
    {
        return m_xPostion;
    }


    enum bcType {none = 0, dirichlet = 1, neumann = 2};

    class boundary_node
    {
    private:
        size_t m_NumNode;
        bcType m_bcType = none;
        double m_BCValue = 0;
    public:
        boundary_node() = default;
        boundary_node(const size_t &NumNode, const bcType &bcType, const double &BCValue);
        ~boundary_node() = default;

        auto GetNode();
        auto GetBCType();
    };

    boundary_node::boundary_node (const size_t &NumNode, const bcType &bcType, const double &BCValue) : m_NumNode(NumNode), m_bcType(bcType),
        m_BCValue(BCValue) {};

    auto boundary_node::GetNode()
    {
        return m_NumNode;
    }

    auto boundary_node::GetBCType()
    {
        return m_bcType;
    }

};