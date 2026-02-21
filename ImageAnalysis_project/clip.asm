; EBX will hold h
; [EBP] will hold s
; [EBP+4] will hold x
; we used fixed point arithemetic with 1.0 represented at 0x10000, in which 
; (1) multiplication results must be divided by 0x10000
; (2) dividends must be multiplied by 0x10000 before the division
; the code implements the clip function defined as
;     CLIP(x, s, h) =
;    0.0             if x < s
;    1.0             if x > h
;    (x-s)/(h-s)     otherwise
; it must be guaranteed that h <> s otherwise there will be a division by zero

section .text
global _start

_start:     
    push    esi                 ; callee save
    push    edi                 ; callee save
                                ; ebx = h
    mov     esi, [ebp]          ; esi = s
    mov     edi, [ebp+4]        ; edi = x
    cmp     edi, esi            ; fall through if x < s
    jae     .l1
    xor     ebx, ebx            ; ebx = 0.0
    jmp     .done 
.l1:
    cmp     edi, ebx            ; fall through if x >= h
    jb      .l2
    mov     ebx, 0x10000        ; ebx = 1.0
    jmp     .done
.l2:
    mov     eax, edi            ; eax = x
    sub     eax, esi            ; eax = x-s
    shl     eax, 16             ; eax = (x-s)*0x10000  (left shift by 16 = multiply by 0x10000)
                                ; valid since x-s fits in 16 bits (s <= x < h <= 0xffff)
    xor     edx, edx            ; clear EDX: div uses EDX:EAX as 64-bit dividend
    sub     ebx, esi            ; ebx = h-s
    div     ebx                 ; eax = (x-s)/(h-s)
    mov     ebx, eax            ; ebx = (x-s)/(h-s)
.done:
    pop     edi                 ; callee restore
    pop     esi                 ; callee restore
    ; Exit (Linux int 0x80 syscall)
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80    