; Fixed point arithemetic using signed integers with 1.0 represented at 0x10000
; (1) multiplication results must be divided by 0x10000
; (2) dividends must be multiplied by 0x10000 before the division
; the code implements the mid-tones function defined as
;    MID(x, m) =
;    0.0             if x = 0.0
;    0.5             if x = m
;    1.0             if x = 1.0
;    (m-1)*x / ((2m-1)*x-m) otherwise
;    in practice the special cases do not need to be tested becuase the formula gives the same results
; EBX will hold m, 0 <= m < 1, in fixed point format
; [EBP] will hold x, 0 <= x < 1, in fixed point format

section .text
global _start

_start:     
    push    esi                 ; callee save
    push    edi                 ; callee save
                                ; ebx = m
    mov     edi, [ebp]          ; edi = x
    ; compute the denominator
    mov     ecx, ebx            ; ecx = m
    shl     ecx, 1              ; ecx = 2*m
    sub     ecx, 0x10000        ; ecx = 2*m - 1.0, recall 1.0 = 0x10000
    mov     eax, edi            ; eax = x
    imul    ecx                 ; edx:eax = (2m-1) * x  (signed: 2m-1 is negative when m < 0.5)
    shrd    eax, edx, 16        ; eax = (2m-1)*x
    sub     eax, ebx            ; eax = (2m-1)*x - m
    mov     ecx, eax            ; ecx = (2m-1)*x - m
    ; compute the numerator
    mov     esi, ebx            ; esi = m
    sub     esi, 0x10000        ; esi = m - 0x10000  (m-1 in fixed point; negative since m < 1)
    mov     eax, edi            ; eax = x
    imul    esi                 ; edx:eax = (m-1)*x
    ; shr/shl by 0x10000 cancel: keep raw product as numerator for the division
    ; perform the division
    idiv    ecx                 ; eax = (m-1)*x / (2m-1)*x - m
    mov     ebx, eax            ; ebx = result
.done:
    pop     edi                 ; callee restore
    pop     esi                 ; callee restore
    ; Exit (Linux int 0x80 syscall)
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80    