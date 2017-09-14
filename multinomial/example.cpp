#include <iostream>
#include "multinomial-layered.cpp"



int main()
{
    double probs[] = {0.1, 0.3, 0.6};
    LayeredMarginal LM(Marginal(probs, 3, 1000));


    double current_threshold = -3.0 ;
    Summator s;
    std::cout << "AAA" << std::endl;
    while(s.get() < 0.9999)
    {
        LM.extend(current_threshold);
        std::cout << "AAA" << s.get() << LM.next() <<  std::endl;
        while(s.get() < 0.9999 and LM.next())
        {
            std::cout << exp(LM.get_lProb()) << "\t" << s.get() << std::endl;
            s.add(exp(LM.get_lProb()));
        }
        current_threshold -= 3.0;
    }
}
