_qa(nam)      = :( abstract $(esc(nam)) )
_qa(nam, sup) = :( abstract $(esc(nam)) <: $(esc(sup)) )

_qs(nam)      = :( immutable $(esc(nam)) ; end )
_qs(nam, sup) = :( immutable $(esc(nam)) <: $(esc(sup)) ; end )
