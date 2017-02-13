_qa(nam)      = :( abstract type $(esc(nam)) end )
_qa(nam, sup) = :( abstract type $(esc(nam)) <: $(esc(sup)) end )

_qs(nam)      = :( struct $(esc(nam)) ; end )
_qs(nam, sup) = :( struct $(esc(nam)) <: $(esc(sup)) ; end )
