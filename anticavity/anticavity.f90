subroutine anticavity(outmask, reflectivity, mask, window_len, refl_thld, nx, ny)

integer(kind=4), dimension(ny, nx), intent(out) :: outmask
integer(kind=4), dimension(ny, nx), intent(in) :: mask
real(kind=8), dimension(ny, nx), intent(in) :: reflectivity
integer(kind=4), intent(in) :: window_len, nx, ny
real(kind=8), intent(in) :: refl_thld

integer(kind=4) :: i, j, k, w, cnt, thrld 
real(kind=8) :: mean_refl

integer(kind=4), dimension(window_len) :: window
real(kind=8), dimension(window_len) :: refl_vec

thrld = window_len / 3
cnt = 0

do i = 1, ny
    do j = 1, nx
        outmask(i, j) = mask(i, j)
    enddo
enddo

do i = 1, ny
    do j = 1, nx - window_len
        window = mask(i, j: j + window_len)
        w = SUM(window)
        if (w == window_len) then
            continue
        elseif (w == 0) then
            continue
        elseif (w > thrld) then
            continue
        endif

        refl_vec = reflectivity(i, j: j + window_len)
        mean_refl = MEAN(refl_vec)

        do k = 0, window_len - 1
            if(window[k + 1] == 1) then
                continue
            endif
            if(ABS(refl_vec[k + 1] - mean_refl) <= refl_thld) then
                outmask(i, j + k) = 1
                cnt = cnt + 1
            endif
        enddo
    enddo
enddo

end subroutine anticavity