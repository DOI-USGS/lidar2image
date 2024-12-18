#ifndef ATK_COMMON_NOT_COPYABLE_H
#define ATK_COMMON_NOT_COPYABLE_H

//To prevent an implicit copy constructor being added to your class that will
//wreak havoc if there is nontrivial state that cannot simply be copied
//(bare pointers are a classical instance), inherit from this class.

//FIXME: Is there a better generally accepted way to do this nowadays?

namespace at
{

class NotCopyable
{
  public:
    NotCopyable() { }
    ~NotCopyable() { }

  private:
    NotCopyable(const NotCopyable &other) { }
};

};

#endif /* ATK_COMMON_NOT_COPYABLE_H */

