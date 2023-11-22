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



    class boundary_node :  public node 
    {
    private:
        int m_bcType = 0;  // usar enum
    public:
        boundary_node() = default;
        boundary_node(const double &x_, const int &bcType_);
        ~boundary_node() = default;
    };

    boundary_node::boundary_node (const double &x_, const int &bcType_) :  node(x_), m_bcType(bcType_) {};

};