Arrays and numbers
==================

While going through the documentation, you may come across the "array" or "number" terms to refer to some function argument or return type.
Arrays and numbers are indeed not proper types. Numbers are any objects satisfying the strengths.typechecking.isnumber() function,
while arrays are any objects satisfying the strengths.typechecking.isarray() function. This also means that no particular type, aside from this restriction,
should be expected. As an example, if a function precises it returns a int, the user can be confident that the returned object will indeed be an int, however,
if the function is said to return a number, the user should not expect any specific type, as it could be a int or a float for instance.
