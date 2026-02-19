; EBX will contain the number of pixels in an image
; [EBP] will contain the address of the histogram (size 0x10000 bytes) prepared by histogram.asm
; the histogram is an array of 32 bit longwords each containing the count of pixels in an image having that respective brighness level
; the pixels in an image each have a brightness value 0...0xffff
; the first longword counts the number of pixels having a value of zero
; last longword in the buffer contains the count of pixels having value 0xffff
; EAX, ECX, EDX are scratch registers that do not need to be preserved
; if other registers are used they must be preserved

; this code computes the median pixel value (we do not distinguish lower and upper medians)

section .text
global _start

_start:
    shr     ebx, 1      ; divide the number of pixels by 2, this is the target we need to reach
    xor     edx, edx    ; edx will count through the number of histogram bins
                        ; xor reg, reg is the idiomatic way to zero a register (shorter and faster than mov reg, 0)
    mov     ecx, [ebp]  ; ecx will be the address of the current histogram bin

.loop:
    sub     ebx, dword [ecx] ; subtract the number of pixels in the current bin from the remaining target
                            ; sub sets CF on unsigned underflow (borrow), so jb is correct here - no cmp needed
    jb      .done           ; edx contains the median pixel value
    inc     edx             ; advance to the next bin
    add     ecx, 4          ; advance to the address of the next bin
    
    cmp     edx, 0x10000    ; test if edx < 0x10000, i.e. there are still bins remaining
    jb      .loop           
    
    or      edx, -1         ; set edx = -1
    
.done                       
    ; Exit (Linux int 0x80 syscall)
    mov     eax, 1              ; sys_exit
    xor     ebx, ebx            ; status 0
    int     0x80
