#ifndef SYSTEMFUNCTOR_HPP  // guard
#define SYSTEMFUNCTOR_HPP 


template <class System, class State, class Time>
class SystemFunctor {
   private:
    System system;
    State parameters;

   public:
    SystemFunctor(System _system, State _parameters)
        : system(_system), parameters(_parameters){};

    void operator()(const State &u, State &dudt, const Time &t) {
        this->system(u, dudt, parameters, t);
    };
};

#endif