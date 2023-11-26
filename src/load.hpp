#pragma once


namespace  FiniteElement
{
    class ExternalLoad
    {
    private:
        size_t NumElem = 0;
        double LoadValue = 0;
    public:
        ExternalLoad() = default;
        ExternalLoad(const size_t &NumElem_, double LoadValue_);
        ~ExternalLoad() = default;


        double GetExternalLoad();
        size_t GetElem();

    };
};


FiniteElement::ExternalLoad::ExternalLoad(const size_t &NumElem_, double LoadValue_): 
                NumElem(NumElem_), LoadValue(LoadValue_){};

double FiniteElement::ExternalLoad::GetExternalLoad()
{
    return LoadValue;
};

size_t FiniteElement::ExternalLoad::GetElem()
{
    return NumElem;
}