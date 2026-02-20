; EBX will contain a number of pixels, n
; [EBP] will contain the address of a buffer of size 0x40000, which is guaranteed to be zero intitialized
; [EBP+4] will contain the address of an image array of n 16-bit words, each containing an unsigned value 0..0xffff, representing the brightness of a pixel
; [EBP+8] will contain the median pixel value of the image array (previously calculated)
; EAX, ECX, EDX are scratch registers that do not need to be preserved
; if other registers are used they must be preserved

; This code converts creates a histogram of the absolute deviaton of pixel brightness values in the buffer from the median
; Each 32-bit longword in the histogram buffer contains the count of the number of pixels having that absolute deviation
; The first longword contains the count of pixels having an absolute deviation of 0x0000 and the
; last longword in the buffer contains the count of pixels having an absolute deviation 0xffff

section .text
global _start

_start:
    push    edi                 ; callee save
    push    esi                 ; callee save (ESI holds image pointer, freeing EDX for CDQ)
    mov     edi, [ebp+8]        ; median pixel value
    mov     esi, [ebp+4]        ; address of the image array
    mov     ecx, [ebp]          ; pointer to the base of the histogram buffer
    test    ebx, ebx            ; exit if no pixels to process
    jz      .done

.loop:
    movzx   eax, word [esi]     ; load a 16-bit word from the image buffer and zero extend to 32 bits
    sub     eax, edi            ; diff = pixel - median; EDX free because ESI holds image pointer
    cdq                         ; EDX = 0x00000000 if EAX >= 0, 0xFFFFFFFF if EAX < 0
    xor     eax, edx            ; if negative: flip all bits
    sub     eax, edx            ; if negative: add 1  ->  EAX = abs(pixel - median)
    inc     dword [ecx + eax*4] ; increment the histogram bin for this absolute deviation
    add     esi, 2              ; move source pointer forward 2 bytes
    dec     ebx                 ; decrement pixel counter
    jnz     .loop               ; continue if so
    
.done:
    ; Exit (Linux int 0x80 syscall)
    pop     esi                 ; callee restore
    pop     edi                 ; callee restore
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80