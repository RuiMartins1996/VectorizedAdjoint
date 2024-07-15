#ifndef BUTCHERTABLE_HPP_INCLUDED
#define BUTCHERTABLE_HPP_INCLUDED

// Boost
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

// TODO: Can use boost::odeint coefficients instead of doing this by hand
// TODO: Could access coefficients directly like so :
// boost::numeric::odeint::rk4_coefficients_a1();

class ButcherTable {
   private:
    std::unique_ptr<std::vector<double>> p_a;
    std::unique_ptr<std::vector<double>> p_b;
    std::unique_ptr<std::vector<double>> p_c;

   public:
    template <typename Stepper, typename State>
    ButcherTable(Stepper stepper,State &u0) {
        typedef typename boost::numeric::odeint::unwrap_reference<
            Stepper>::type::stepper_category stepper_category;

        initialize(stepper, u0, stepper_category());
    };

    template <typename Stepper, typename State>
    void initialize(Stepper stepper,State &u0,
                    boost::numeric::odeint::explicit_controlled_stepper_tag) {
        std::cout << "To create Driver, we need to supply and error stepper or a regular stepper, not a controlled stepper!" << std::endl;
    }

    template <typename Stepper, typename State>
    void initialize(Stepper stepper, State &u0,
                    boost::numeric::odeint::stepper_tag) {
        if ((typeid(stepper) == typeid(boost::numeric::odeint::euler<State>))) {
            //std::cout << "Euler\n";

            int stages = 1;
            std::unique_ptr<std::vector<double>> _a =
                std::make_unique<std::vector<double>>(stages * stages);
            std::unique_ptr<std::vector<double>> _b =
                std::make_unique<std::vector<double>>(stages);

            std::unique_ptr<std::vector<double>> _c =
                std::make_unique<std::vector<double>>(stages);

            for (int i = 0; i < stages * stages; i++) (*_a)[i] = 0.0;
            (*_b)[0] = 1.0;
            (*_c)[0] = 0.0;

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        }

        if (typeid(stepper) ==
            typeid(boost::numeric::odeint::runge_kutta4_classic<State>)) {
            //std::cout << "Runge-Kutta 4\n";

            int stages = 4;

            std::unique_ptr<std::vector<double>> _a =
                std::make_unique<std::vector<double>>(stages * stages);
            std::unique_ptr<std::vector<double>> _b =
                std::make_unique<std::vector<double>>(stages);

            std::unique_ptr<std::vector<double>> _c =
                std::make_unique<std::vector<double>>(stages);

            for (int i = 0; i < stages * stages; i++) (*_a)[i] = 0.0;

            auto id = [stages](int i, int j) { return i * stages + j; };

            (*_a)[id(1, 0)] = static_cast<double>(1) / static_cast<double>(2);
            (*_a)[id(2, 1)] = static_cast<double>(1) / static_cast<double>(2);
            (*_a)[id(3, 2)] = static_cast<double>(1) / static_cast<double>(2);

            (*_b)[0] = static_cast<double>(1) / static_cast<double>(6);
            (*_b)[1] = static_cast<double>(1) / static_cast<double>(3);
            (*_b)[2] = static_cast<double>(1) / static_cast<double>(3);
            (*_b)[3] = static_cast<double>(1) / static_cast<double>(6);

            (*_c)[0] = static_cast<double>(0);
            (*_c)[1] = static_cast<double>(1) / static_cast<double>(2);
            (*_c)[2] = static_cast<double>(1) / static_cast<double>(2);
            (*_c)[3] = static_cast<double>(1);

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        };

        if (typeid(stepper) ==
            typeid(boost::numeric::odeint::explicit_generic_rk<3,3,State>)) {
                std::cout<< "Hello\n";
            }
    };

    template <typename Stepper, typename State>
    void initialize(Stepper stepper, State &u0,
                    boost::numeric::odeint::error_stepper_tag) {
        if (typeid(stepper) ==
            typeid(boost::numeric::odeint::runge_kutta_cash_karp54<State>)) {
            //std::cout << "Runge-Kutta CashKarp(45) 4\n";

            int stages = 6;
            std::unique_ptr<std::vector<double>> _a =
                std::make_unique<std::vector<double>>(stages * stages);
            std::unique_ptr<std::vector<double>> _b =
                std::make_unique<std::vector<double>>(stages);

            std::unique_ptr<std::vector<double>> _c =
                std::make_unique<std::vector<double>>(stages);

            for (int i = 0; i < stages * stages; i++) (*_a)[i] = 0.0;

            auto id = [stages](int i, int j) { return i * stages + j; };

            (*_a)[id(1, 0)] = static_cast<double>(1) / static_cast<double>(5);
            (*_a)[id(2, 0)] = static_cast<double>(3) / static_cast<double>(40);
            (*_a)[id(2, 1)] = static_cast<double>(9) / static_cast<double>(40);
            (*_a)[id(3, 0)] = static_cast<double>(3) / static_cast<double>(10);
            (*_a)[id(3, 1)] = static_cast<double>(-9) / static_cast<double>(10);
            (*_a)[id(3, 2)] = static_cast<double>(6) / static_cast<double>(5);
            (*_a)[id(4, 0)] =
                static_cast<double>(-11) / static_cast<double>(54);
            (*_a)[id(4, 1)] = static_cast<double>(5) / static_cast<double>(2);
            (*_a)[id(4, 2)] =
                static_cast<double>(-70) / static_cast<double>(27);
            (*_a)[id(4, 3)] = static_cast<double>(35) / static_cast<double>(27);
            (*_a)[id(5, 0)] =
                static_cast<double>(1631) / static_cast<double>(55296);
            (*_a)[id(5, 1)] =
                static_cast<double>(175) / static_cast<double>(512);
            (*_a)[id(5, 2)] =
                static_cast<double>(575) / static_cast<double>(13824);
            (*_a)[id(5, 3)] =
                static_cast<double>(44275) / static_cast<double>(110592);
            (*_a)[id(5, 4)] =
                static_cast<double>(253) / static_cast<double>(4096);

            (*_b)[0] = static_cast<double>(37) / static_cast<double>(378);
            (*_b)[1] = static_cast<double>(0);
            (*_b)[2] = static_cast<double>(250) / static_cast<double>(621);
            (*_b)[3] = static_cast<double>(125) / static_cast<double>(594);
            (*_b)[4] = static_cast<double>(0);
            (*_b)[5] = static_cast<double>(512) / static_cast<double>(1771);

            (*_c)[0] = static_cast<double>(0);
            (*_c)[1] = static_cast<double>(1) / static_cast<double>(5);
            (*_c)[2] = static_cast<double>(3) / static_cast<double>(10);
            (*_c)[3] = static_cast<double>(3) / static_cast<double>(5);
            (*_c)[4] = static_cast<double>(1);
            (*_c)[5] = static_cast<double>(7) / static_cast<double>(8);

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        };

        if (typeid(stepper) ==
            typeid(boost::numeric::odeint::runge_kutta_fehlberg78<State>)) {
            //std::cout << "Runge-Kutta Fehlberg(78) 4\n";

            int stages = 13;
            std::unique_ptr<std::vector<double>> _a =
                std::make_unique<std::vector<double>>(stages * stages);
            std::unique_ptr<std::vector<double>> _b =
                std::make_unique<std::vector<double>>(stages);

            std::unique_ptr<std::vector<double>> _c =
                std::make_unique<std::vector<double>>(stages);

            for (int i = 0; i < stages * stages; i++) (*_a)[i] = 0.0;

            auto id = [stages](int i, int j) { return i * stages + j; };

            (*_a)[id(1, 0)] = static_cast<double>(2) / static_cast<double>(27);

            (*_a)[id(2, 0)] = static_cast<double>(1) / static_cast<double>(36);
            (*_a)[id(2, 1)] = static_cast<double>(1) / static_cast<double>(12);

            (*_a)[id(3, 0)] = static_cast<double>(1) / static_cast<double>(24);
            (*_a)[id(3, 1)] = static_cast<double>(0);
            (*_a)[id(3, 2)] = static_cast<double>(1) / static_cast<double>(8);

            (*_a)[id(4, 0)] = static_cast<double>(5) / static_cast<double>(12);
            (*_a)[id(4, 1)] = static_cast<double>(0);
            (*_a)[id(4, 2)] =
                static_cast<double>(-25) / static_cast<double>(16);
            (*_a)[id(4, 3)] = static_cast<double>(25) / static_cast<double>(16);

            (*_a)[id(5, 0)] = static_cast<double>(1) / static_cast<double>(20);
            (*_a)[id(5, 1)] = static_cast<double>(0);
            (*_a)[id(5, 2)] = static_cast<double>(0);
            (*_a)[id(5, 3)] = static_cast<double>(1) / static_cast<double>(4);
            (*_a)[id(5, 4)] = static_cast<double>(1) / static_cast<double>(5);

            (*_a)[id(6, 0)] =
                static_cast<double>(-25) / static_cast<double>(108);
            (*_a)[id(6, 1)] = static_cast<double>(0);
            (*_a)[id(6, 2)] = static_cast<double>(0);
            (*_a)[id(6, 3)] =
                static_cast<double>(125) / static_cast<double>(108);
            (*_a)[id(6, 4)] =
                static_cast<double>(-65) / static_cast<double>(27);
            (*_a)[id(6, 5)] =
                static_cast<double>(125) / static_cast<double>(54);

            (*_a)[id(7, 0)] =
                static_cast<double>(31) / static_cast<double>(300);
            (*_a)[id(7, 1)] = static_cast<double>(0);
            (*_a)[id(7, 2)] = static_cast<double>(0);
            (*_a)[id(7, 3)] = static_cast<double>(0);
            (*_a)[id(7, 4)] =
                static_cast<double>(61) / static_cast<double>(225);
            (*_a)[id(7, 5)] = static_cast<double>(-2) / static_cast<double>(9);
            (*_a)[id(7, 6)] =
                static_cast<double>(13) / static_cast<double>(900);

            (*_a)[id(8, 0)] = static_cast<double>(2);
            (*_a)[id(8, 1)] = static_cast<double>(0);
            (*_a)[id(8, 2)] = static_cast<double>(0);
            (*_a)[id(8, 3)] = static_cast<double>(-53) / static_cast<double>(6);
            (*_a)[id(8, 4)] =
                static_cast<double>(704) / static_cast<double>(45);
            (*_a)[id(8, 5)] =
                static_cast<double>(-107) / static_cast<double>(9);
            (*_a)[id(8, 6)] = static_cast<double>(67) / static_cast<double>(90);
            (*_a)[id(8, 7)] = static_cast<double>(3);

            (*_a)[id(9, 0)] =
                static_cast<double>(-91) / static_cast<double>(108);
            (*_a)[id(9, 1)] = static_cast<double>(0);
            (*_a)[id(9, 2)] = static_cast<double>(0);
            (*_a)[id(9, 3)] =
                static_cast<double>(23) / static_cast<double>(108);
            (*_a)[id(9, 4)] =
                static_cast<double>(-976) / static_cast<double>(135);
            (*_a)[id(9, 5)] =
                static_cast<double>(311) / static_cast<double>(54);
            (*_a)[id(9, 6)] =
                static_cast<double>(-19) / static_cast<double>(60);
            (*_a)[id(9, 7)] = static_cast<double>(17) / static_cast<double>(6);
            (*_a)[id(9, 8)] = static_cast<double>(-1) / static_cast<double>(12);

            (*_a)[id(10, 0)] =
                static_cast<double>(2383) / static_cast<double>(4100);
            (*_a)[id(10, 1)] = static_cast<double>(0);
            (*_a)[id(10, 2)] = static_cast<double>(0);
            (*_a)[id(10, 3)] =
                static_cast<double>(-341) / static_cast<double>(164);
            (*_a)[id(10, 4)] =
                static_cast<double>(4496) / static_cast<double>(1025);
            (*_a)[id(10, 5)] =
                static_cast<double>(-301) / static_cast<double>(82);
            (*_a)[id(10, 6)] =
                static_cast<double>(2133) / static_cast<double>(4100);
            (*_a)[id(10, 7)] =
                static_cast<double>(45) / static_cast<double>(82);
            (*_a)[id(10, 8)] =
                static_cast<double>(45) / static_cast<double>(164);
            (*_a)[id(10, 9)] =
                static_cast<double>(18) / static_cast<double>(41);

            (*_a)[id(11, 0)] =
                static_cast<double>(3) / static_cast<double>(205);
            (*_a)[id(11, 1)] = static_cast<double>(0);
            (*_a)[id(11, 2)] = static_cast<double>(0);
            (*_a)[id(11, 3)] = static_cast<double>(0);
            (*_a)[id(11, 4)] = static_cast<double>(0);
            (*_a)[id(11, 5)] =
                static_cast<double>(-6) / static_cast<double>(41);
            (*_a)[id(11, 6)] =
                static_cast<double>(-3) / static_cast<double>(205);
            (*_a)[id(11, 7)] =
                static_cast<double>(-3) / static_cast<double>(41);
            (*_a)[id(11, 8)] = static_cast<double>(3) / static_cast<double>(41);
            (*_a)[id(11, 9)] = static_cast<double>(6) / static_cast<double>(41);
            (*_a)[id(11, 10)] = static_cast<double>(0);

            (*_a)[id(12, 0)] =
                static_cast<double>(-1777) / static_cast<double>(4100);
            (*_a)[id(12, 1)] = static_cast<double>(0);
            (*_a)[id(12, 2)] = static_cast<double>(0);
            (*_a)[id(12, 3)] =
                static_cast<double>(-341) / static_cast<double>(164);
            (*_a)[id(12, 4)] =
                static_cast<double>(4496) / static_cast<double>(1025);
            (*_a)[id(12, 5)] =
                static_cast<double>(-289) / static_cast<double>(82);
            (*_a)[id(12, 6)] =
                static_cast<double>(2193) / static_cast<double>(4100);
            (*_a)[id(12, 7)] =
                static_cast<double>(51) / static_cast<double>(82);
            (*_a)[id(12, 8)] =
                static_cast<double>(33) / static_cast<double>(164);
            (*_a)[id(12, 9)] =
                static_cast<double>(12) / static_cast<double>(41);
            (*_a)[id(12, 10)] = static_cast<double>(0);
            (*_a)[id(12, 11)] = static_cast<double>(1);

            (*_b)[0] = static_cast<double>(0);
            (*_b)[1] = static_cast<double>(0);
            (*_b)[2] = static_cast<double>(0);
            (*_b)[3] = static_cast<double>(0);
            (*_b)[4] = static_cast<double>(0);
            (*_b)[5] = static_cast<double>(34) / static_cast<double>(105);
            (*_b)[6] = static_cast<double>(9) / static_cast<double>(35);
            (*_b)[7] = static_cast<double>(9) / static_cast<double>(35);
            (*_b)[8] = static_cast<double>(9) / static_cast<double>(280);
            (*_b)[9] = static_cast<double>(9) / static_cast<double>(280);
            (*_b)[10] = static_cast<double>(0);
            (*_b)[11] = static_cast<double>(41) / static_cast<double>(840);
            (*_b)[12] = static_cast<double>(41) / static_cast<double>(840);

            (*_c)[0] = static_cast<double>(0);
            (*_c)[1] = static_cast<double>(2) / static_cast<double>(27);
            (*_c)[2] = static_cast<double>(1) / static_cast<double>(9);
            (*_c)[3] = static_cast<double>(1) / static_cast<double>(6);
            (*_c)[4] = static_cast<double>(5) / static_cast<double>(12);
            (*_c)[5] = static_cast<double>(1) / static_cast<double>(2);
            (*_c)[6] = static_cast<double>(5) / static_cast<double>(6);
            (*_c)[7] = static_cast<double>(1) / static_cast<double>(6);
            (*_c)[8] = static_cast<double>(2) / static_cast<double>(3);
            (*_c)[9] = static_cast<double>(1) / static_cast<double>(3);
            (*_c)[10] = static_cast<double>(1);
            (*_c)[11] = static_cast<double>(0);
            (*_c)[12] = static_cast<double>(1);

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        };
    };

   public:
    int GetStages() { return (*p_b).size(); };
    double a(int i, int j) {
        int s = GetStages();
        return (*p_a)[s * i + j];
    };

    double b(int i) { return (*p_b)[i]; };

    double c(int i) { return (*p_c)[i]; }
};

#endif
