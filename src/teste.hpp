#include <vector>

class BaseType
{
public:
    virtual BaseType*  myfunction() = 0;
    virtual ~BaseType() {}
};

class IntType : public BaseType
{
    int X;
    BaseType*  myfunction();
};

class BoolType  : public BaseType
{
    bool b;
    BaseType*  myfunction();
};

class CharType : public BaseType
{
    char c;
    BaseType*  myfunction();
};

BaseType*  myfunction(BaseType* b)
{

    return b->myfunction();
}


BaseType a;
a->myfunction();

std::vector<BaseType> name;