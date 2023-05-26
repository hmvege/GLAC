#include "action.h"
#include "config/parameters.h"

Action::Action()
{
    m_N = Parameters::getN();
    m_position = std::vector<int>(4,0);
}

Action::~Action()
{
}
