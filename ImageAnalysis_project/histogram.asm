; EBX will contain a number of pixels, n.  It is guaranteed that n is greater than zero
; [EBP] will contain the address of a buffer of size 0x40000, which is guaranteed to be zero intitialized
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
    test    ebx, ebx        ; check if pixel count is zero
    jz      .done           ; exit if no pixels to process
.loop:
    movzx   eax, word [edx]     ; load 16-bit word and zero extend to 32 bits
    inc     dword [ecx + eax*4] ; increment the longword at address = ax + ecx
    
    ; Advance pointers and decrement counter
    add     edx, 2              ; move source pointer forward 2 bytes
    dec     ebx                 ; decrement pixel count by 1
    jnz     .loop               ; continue if ebx > 0

.done:
    ; Exit (Linux int 0x80 syscall)
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80