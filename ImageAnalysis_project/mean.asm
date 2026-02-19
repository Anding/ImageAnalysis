; EBX will contain the number of pixels in an image
; [EBP] will contain the address of the histogram (size 0x40000 bytes) prepared by histogram.asm
; the histogram is an array of 32 bit longwords each containing the count of pixels in an image having that respective brighness level
; the pixels in an image each have a brightness value 0...0xffff
; the first longword counts the number of pixels having a value of zero
; last longword in the buffer contains the count of pixels having value 0xffff
; EAX, ECX, EDX are scratch registers that do not need to be preserved
; if other registers are used they must be preserved

; this code computes the mean pixel value

section .text
global _start

_start:
    push    esi             ; callee save
    push    edi             ; callee save
    push    ebx             ; save the number of pixels on the stack for the final division operation
    xor     esi, esi        ; esi will be the hi 32 bits of a 64-bit accumulator
    xor     edi, edi        ; edi will be the lo 32 bits of a 64-bit accumulator
    xor     ecx, ecx        ; ecx will be the number of the current histogram bin 0..0xffff
    mov     ebx, [ebp] ; ebx will be the address of the current histogram bin

.loop:
    mov     eax, ecx        ; eax = the number of the current histogram bin
    mul     dword [ebx]     ; edx:eax = eax * [ebx], the total pixel intensity in this bin
    add     edi, eax        ; lo 32 bits of a 64-bit addition
    adc     esi, edx        ; hi 32 bits of a 64-bit addition
    inc     ecx             ; advance to the next bin
    add     ebx, 4          ; advance to the address of the next bin
    cmp     ecx, 0x10000   ; test if ecx < 0x10000, i.e. there are still bins remaining
    jb      .loop           
    
    mov     edx, esi        ; move the accumulated pixel intensity to edx:eax
    mov     eax, edi
    pop     ebx             ; reload ebx with the number of pixels
    div     ebx             ; eax now contains the mean pixel value
    pop     edi             ; callee restore
    pop     esi             ; callee restore

.done                       
    ; Exit (Linux int 0x80 syscall)
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80
