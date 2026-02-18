; EBX will contain a number of bytes, n
; [EBP] will contain the address of a buffer of size 0x4000, which is guaranteed to be zero intitialized
; [EBP+4] will contain the address of an array of n 16-bit words, each containing an unsigned value 0..0xffff, representing the brightness of a pixel
; EAX, ECX, EDX are scratch registers that do not need to be preserved
; if other registers are used they must be preserved

; This code converts creates a histogram of pixel brightness values in the buffer.  Each 32-bit longword in the buffer contains the
; count of the number of pixels having that brightness value.  The first longword contains the count of pixels having a value of 0x0000 and the
; last longword in the buffer contains the count of pixels having value 0xffff

section .text
global _start

_start:
    ; Load pointers: EDX = array of image pixels (16-bit words), ECX = histogram buffer (32-bit longwords). EBX = number of bytes
    mov     edx, [ebp+4]    ; source pointer
    mov     ecx, [ebp]      ; destination pointer

    test    ebx, ebx        ; check if byte count is zero
    jz      .done

.loop:
    cmp     ebx, 2          ; check if at least 2 bytes remain
    jb      .done
    
    ; Process one 16-bit word
    mov     ax, word [edx]      ; load 16-bit word
    inc     dword [ecx + eax*4] ; increment the longword at address = ax + ecx
    
    ; Advance pointers and decrement counter
    add     edx, 2              ; move source pointer forward 2 bytes
    sub     ebx, 2              ; decrement byte counter by 2
    jmp     .loop               ; continue loop

.done:
    ; Exit (Linux int 0x80 syscall)
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80