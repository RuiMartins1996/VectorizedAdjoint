#ifndef STATE_STORAFE_NEW_HPP  // guard
#define STATE_STORAFE_NEW_HPP


struct StateStorage {
   private:
    std::vector<double> time;
    std::vector<std::vector<double>> states;

   public:
    StateStorage() = default;
    ~StateStorage() {
        //std::cout << "Calling destructor for StateStorage" << std::endl;
    };
    void Clear() {
        time.clear();
        states.clear();
    }
    int GetT() const { return time.size(); }
    double GetTime(const int n) const { return time[n]; };
    double GetDt(const int n) const { return time[n + 1] - time[n]; }
    template <typename Time>
    void PushBackTime(const Time _time) {
        time.push_back(static_cast<double>(_time));
    }
    template <typename State>
    void PushBackState(State state) {
        std::vector<double> state_push_back(state.size());
        for (int i = 0; i < state.size(); i++) state_push_back[i] = state[i];
        states.push_back(state_push_back);
    }

    template <class State>
    void GetState(State &u, int n) const {
        for (int i = 0; i < u.size(); i++) u[i] = states[n][i];
    }
};

#endif