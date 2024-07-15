#pragma once

#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include<cstddef>
#include<iostream>
#include<cstring>
#include<stdexcept>
#include<memory>
#include<list>
#include<unordered_map>
#include<typeinfo>
#include<typeindex>
#include <type_traits>

/*
 * How to convert normal data class into AADC compatible version:
 *  1) Derive data class from aadc::VReg as "public virtual aadc::VReg"
 *  2) Add object registration into data class constructor: AADC_VREG_OBJ_ADD;
 *  3) Add virtual table switch function using macro AADC_VREG_SWITCH;
 *
 * 	Example to convert normal class A into AADC version aadc::A:
 *
 * 	AADC_VREG_LINK(aadc::A, A);							// register link between two types: AADC and plain version
 * 	auto *a = new A(1);								// instance of plain class
 * 	auto aa = reinterpret_cast<aadc::A*>(a);        // cast plain object into AADC version
 *
 *  aadc::VReg::SwitchAll();                        // switch virtual tables from normal into aadc version
 *
 * TODO: make direction of the switch explicit instead of implicit
 */

namespace aadc {

    class VReg {


        // record active objects here
        static std::list<VReg*>& actives() { static std::list<VReg*> actives; return actives; }
        // size = obj_size: helper array: indicate if vtable at given object offset has been recorded
        static std::vector<bool>& obj_start_buffer() { static std::vector<bool> obj_start_buffer; return obj_start_buffer; }; // used during obj LINKING
        // size = obj_size: temp buffer to hold first object
        static std::vector<uint8_t>& obj_mem() { static std::vector<uint8_t> obj_mem; return obj_mem; }
        // address of object head : &obj_mem
        static void*& obj_start() { static void *obj_start; return obj_start; }
        // indicate if we are counting vtables for the first type
        static bool& is_linking() { static bool is_linking; return is_linking; }
        // link type's vtable to type index of the linked type, type index is given by std::type_index(typeid(T1))
        typedef std::unordered_map<std::type_index, std::vector<std::pair<int, void*> > > Map;

        static Map& obj_vtables() { static Map obj_vtables; return obj_vtables; }
        // number of vtables for the top object
        static int& num_vtables() { static int num_vtables; return num_vtables; }

    public:
        VReg() : pos_iter(actives().insert(actives().end(), this))
        {
            if(debug()) std::cout << "VReg() " << typeid(*this).name() << std::endl;
        }

        VReg(const VReg& other) : pos_iter(actives().insert(actives().end(), this))
        {}

        virtual ~VReg() {
            actives().erase(pos_iter);
            if(debug()) std::cout << "~VReg() " << typeid(*this).name() << std::endl;
        }

        // need to restore clean state for unit tests
        static void Clean() {
            actives().clear();
            obj_start_buffer().clear();
            obj_mem().clear();
            obj_vtables().clear();
            obj_start() = nullptr;
            num_vtables() = 0;
        }

        // must be defined in the derived client class
        virtual size_t switchVT() = 0;

        template<class T>
        void debugInfo(const T* t) {
            auto name = typeid(T).name();
            auto size = sizeof(T);
            if (debug())
                std::cout << "add obj=" << t << " type=" << name << " actives=" << actives().size() << " vtables=" << num_vtables()  << " objsize=" << size << std::endl;
        }

        // 1. count vtables
        // 2. update obj_start_buffer
        template<class T>
        void AddObj(const T *t) {
            if (!is_linking())                    return;
            if (!std::is_polymorphic<T>::value) return;

			uint64_t offset((uint64_t)t - (uint64_t)obj_start());
            if (!obj_start_buffer()[offset]) ++num_vtables();
            obj_start_buffer()[offset] = true;
        }

        // T1 is an active type & T2 is a plain type
        template<class T1, class T2>
        static void AddLink() {
            static_assert(sizeof(T1) == sizeof(T2), "AADC Objects Linking. Objects are not compatible.");

            VReg::LinkStart(sizeof(T1));

            T1* t1 = new (&(VReg::obj_mem().front())) T1();
            VReg::is_linking() = false;
            T2 t2;
            VReg::LinkVTables(t1, &t2);
            t1->~T1();
        }

        template<class T1, class T2>
        static void LinkVTables(const T1* t1, const T2* t2) {
            std::vector<std::pair<int, void*> > vt1(num_vtables()), vt2(num_vtables());
            for (int i = 0, vti = 0; i < obj_start_buffer().size(); ++i) {
                if (obj_start_buffer()[i]) {
                    vt1[vti] = std::pair<int, void*>(i, *(void**)((uint64_t)t1+i));
                    vt2[vti] = std::pair<int, void*>(i, *(void**)((uint64_t)t2+i));
                    ++vti;
                }
            }
            obj_vtables()[std::type_index(typeid(T1))] = vt2;
            obj_vtables()[std::type_index(typeid(T2))] = vt1;
            is_linking() = false;
        }

        // called from client code
        static void SwitchAll() {
            for(auto & active : actives()) {
                active->switchVT();
            }
        }

        // called from derived class
        template<class T1>
        static size_t SwitchVT(T1* t) {
            std::vector<std::pair<int, void*> >& vt(obj_vtables()[std::type_index(typeid(T1))]); // TODO: check if exists

            for (int vti = 0; vti < vt.size(); ++vti) {
                *(void**)((uint64_t)t + vt[vti].first) = vt[vti].second;
            }
            return vt.size();
        }

        // these are for diagnostics and unit tests only
        static std::list<VReg*>& GetActives()                    { return actives(); }
        static int               GetVTCount()                    { return num_vtables(); }
        static bool& debug() { static bool debug; return debug; }

    private:
        static void LinkStart(int obj_size) {
            obj_mem().resize(obj_size);
            obj_start_buffer().resize(obj_size);
            std::fill(obj_start_buffer().begin(), obj_start_buffer().end(), false);
            obj_start() = &(obj_mem().front());
            is_linking() = true;
            num_vtables() = 0;
        }



        const std::list<VReg*>::iterator pos_iter;
    };
}

// objects must be derived from aadc::VReg
#define AADC_VREG_BASE public virtual aadc::VReg

// must be placed inside market class to enable switch functinality
#define	AADC_VREG_SWITCH virtual size_t switchVT() { return aadc::VReg::SwitchVT(this); }

// object registration macro : has to go into default object constructor
#define	AADC_VREG_OBJ_ADD	{ aadc::VReg::AddObj(this); aadc::VReg::debugInfo(this); }

// object de-registration macro : not needed now, empty
#define	AADC_VREG_OBJ_DEL	; /* aadc::VReg::DelObj(this) */

// register link between two object types: one AADC and one non-AADC
// we don't need to keep instance of the object: we only need to create instance, collect type information and let object go
#define	AADC_VREG_LINK(T1, T2)	aadc::VReg::AddLink<T1, T2>();



