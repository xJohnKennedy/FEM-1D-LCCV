#include <vector>
#include <initializer_list>

#pragma once

namespace FiniteElement {
    
    namespace  Material {

        class Material 
            {
            public:
                
                Material() = default;
                ~Material() = default;
                
                virtual double GetDensity() { return 0; };
                virtual double GetDensity(const double &x_ = 0) { return 0; };
                virtual std::vector<double> GetConstitutiveMatrix(double ParametricPostition){return std::vector<double>{};};
            };


        class Linear : public Material
        {
        protected:

        public:

            Linear() = default;
            ~Linear() = default;
        };

        class Homogeneous : public Linear 
        {
            private:
                double m_density = 0;
                double m_modulus = 0;
            public:

                Homogeneous() = default;
                ~Homogeneous() = default;
                Homogeneous (const double &density_, const double &modulus_);

                double GetDensity();
                double GetDensity(const double &x_);
                std::vector<double> GetConstitutiveMatrix(double ParametricPostition);

        };

        Homogeneous::Homogeneous (const double &density_, const double &modulus_):
            m_density(density_), m_modulus(modulus_) {
            };

        
        double Homogeneous::GetDensity ()
        {
            return m_density;
        };
        
        double Homogeneous::GetDensity (const double &x_ = 0)
        {
            return m_density;
        };

        std::vector<double> Homogeneous::GetConstitutiveMatrix(double ParametricPostition)
        {
            return std::vector<double> {m_modulus};
        };
    };
};