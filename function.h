#ifndef FUNCTION_H
#define FUNCTION_H

class FUNC
{
    //calculate md parameters with temperature, density and natom
    public:
        void mdparameter();

    //convert between different units
    public:
        void unitcovert();
    private:
        void convertT();
    
};
#endif
