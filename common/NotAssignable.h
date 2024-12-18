#ifndef ATK_COMMON_NOT_ASSIGNABLE_H
#define ATK_COMMON_NOT_ASSIGNABLE_H

//To prevent an implicit assignment operator being added to your class that will
//wreak havoc if there is nontrivial state that cannot simply be copied
//(bare pointers are a classical instance), inherit from this class.

//FIXME: Is there a better generally accepted way to do this nowadays?

namespace at
{

class NotAssignable
{
  public:
    NotAssignable() { }
    ~NotAssignable() { }

  private:
    void operator=(const NotAssignable &other) { }
};

};

#endif /* ATK_COMMON_NOT_ASSIGNABLE_H */

