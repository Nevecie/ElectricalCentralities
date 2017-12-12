USE MPI                                                     !mpi module
USE ism                                                     !итерационные решатели СЛАУ (gmres); использует модуль blas 

IMPLICIT NONE                                               !объявлять все

INTEGER:: n, nonzero                                        !размер матрицы и число ненулевых элементов
character(len=25):: fvcr                                    !файл данных, должен содержать 3 строки чисел, см далее
character(len=25):: fres                                    !файл для результата
character(len=66):: fgraph                                  !файл с описанием графа для выдачи ответа в удобной форме
character(len=9), dimension(2):: str                        !вспомогательная текстовая строка для меток графа (нужна, если все ок)
integer:: error                                             !код ошибки
double precision, dimension(:), allocatable:: a             !число ненулевых элементов матрицы или ее верхнего треугольника, слева направо, строка за строкой. Диагональные элементы надо указать, даже если они нули.
double precision, dimension(:), allocatable:: w, CF, CF0    !сохраним в w матрицу смежности, а в CF будут искомые центральности. CF0 соберет их со всех процессоров.
double precision, dimension(:), allocatable:: b, x, bb      !правая часть и будущее решение, а также вспомогательный массив - аккумулятор, буфер и т.п.
integer, dimension(:), allocatable:: ia                     !позиции первых ненулевых элементов строк, один элемент фиктивный. Если элементов нет, берется с предыдущей строки
integer, dimension(:), allocatable:: ja                     !индексы столбцов ненулевых элементов
integer:: iwk,lfil,i,j, smaxres                             !вспомогательные переменные
integer:: Krylov = 7, maxiter = 100                         !параметры для решателя GMRES
DOUBLE PRECISION:: epsss = 1.D-12,droptol                   !параметры для решателя GMRES
double precision, dimension(:,:), allocatable:: vv          !массивы для решателя GMRES
double precision, dimension(:), allocatable:: alu, wa       !массивы для решателя GMRES
integer, dimension(:), allocatable:: ju,jlu,jw              !массивы для решателя GMRES
double precision, dimension(:), allocatable:: D             !диагональ матрицы
double precision:: delta, resid, maxres=-1.,maxres0=-1      !параметр дельта из постановки задачи, максимальная ошибка в решении.
logical:: verbose = .true., check=.true.                    !флаги вывода информации - не для больших задач! и проверки решения системы.
integer:: s,s1,s2,each                                      !индексы вершин - рабочий счетчик и пределы диапазона для процесса, сколько вершин на процесс
integer:: mpierr, comm=MPI_COMM_WORLD, ncpus, myid          !Код ошибки, коммуникатор, число процессов и ранг текущего процесса

NAMELIST /SIZES/ n, nonzero, delta, verbose, check, krylov,epsss, & 
        & lfil,droptol,maxiter,iwk, fgraph, fvcr, fres      !namelist: чтобы считать все значения из файла, см. gmres.nml

call mpi_init(mpierr)                      !Запуск MPI 
call mpi_comm_size(comm,ncpus,mpierr)      !Число процессов
call mpi_comm_rank(comm,myid,mpierr)       !Ранк этого процесса
open(42,file='gmres.nml')                  !открыть файл параметров
read(42,SIZES)                             !читать параметры
close(42)                                  !больше не нужен
iwk = iwk*nonzero                          !константа для решателя

!создать массивы данного размера начиная с 1
allocate(a(nonzero),w(nonzero))            !элементы матриц: рабочей и исходной
allocate(CF(n)); CF = 0.                   !центральности (не полные, а по вершинам данного процесса)
allocate(ja(nonzero))                      !индексы столбцов ненулевых элементов матрицы
allocate(ia(n+1))                          !номера первых ненулевых в каждой строке плюс один фиктивный 
allocate(b(n),bb(n),x(n))                  !правая часть, решение, буфер
allocate(vv(n,Krylov+1))                   !для решателя
allocate(alu(iwk))                         !для решателя и предобуславливателя
allocate(ju(n), wa(n+1), jw(2*n))          !--
allocate(jlu(iwk))                         !--
allocate(D(n))                             !диагональ диагональной матрицы
open(1, file=fvcr)                         !открыть файл; он готовится из 3-колоночного описания графа утилитой *
read(1,*) a                                !читать элементы матрицы
read(1,*) ja                               !читать индексы столбцов
read(1,*) ia                               !читать позиции для различения строк
w=a                                        !сохраним матрицу смежности
close(1)                                   !файл больше не нужен

!готовим матрицу: D - W + delta I, D диагональная, сумма по строке W, I единичная, delta - параметр
do i =1,n
        D(i) = sum( a(ia(i):ia(i+1)-1) )   !сумма по строке
enddo
j=0                                        !счетчик диагональных элементов
do i=1,nonzero                             !преобразовать матрицу: W -> D-W+deltaI. Цикл по списку ненулевых элементов.
        a(i)=-a(i)                         !внедиагональные элементы просто меняют знак
        if(ia(ja(i)).le.i .and. ia(ja(i)+1).gt.i) then !диагональ: элемент принадлежит строке с тем же индексом что и столбец этого элемента
                if((a(i).ne.0.)) print*, 'warning: diagonal nonzero:', i,a(i),ia(i:i+1),ja(i) !диагональ должна быть нулевая, но перечислена вместе с ненулевыми элементами
                j=j+1                      !следующий диагональный элемент
                a(i) = D(j)+a(i)+delta*1.  !фактически a(i)=D(j)+delta т.к. a(i) должны быть нулями, это диагональные элементы, и плюс дельта на единичную матрицу - дельта к диагонали
        endif
enddo
call ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,wa,jw,error)!преобразование матрицы с целью улучшения сходимости метода
select case(error)                         !анализируем код ошибки
case(0)                                    !все нормально
        continue                           !продолжаем работу
case(1:)  
        print*, 'Preconditioning ERROR: zero pivot encountered at step number', error
        call bye
case(-1)
        print*, 'Preconditioning ERROR: input matrix may be wrong. (The elimination process has generated a row in L or U whose length is > n.)'
        call bye
case(-2)
        print*, 'Preconditioning ERROR: The matrix L overflows the array al. Increase its size.'
        call bye
case(-3)
        print*, 'Preconditioning ERROR: The matrix U overflows the array alu. Increase its size iwk.'
        call bye
case(-4)
        print*, 'Preconditioning ERROR: Illegal value for lfil:', lfil
        call bye
case(-5)
        print*, 'Preconditioning ERROR: Zero row encountered.'
        call bye
case default
        print*, 'Unknown Preconditioning ERROR:', error
        call bye
end select

if(verbose) open(65,file='solutions.dat')  !откроем файл для ответа, если включен вывод
each = ceiling(real(n)/real(ncpus))        !работа каждого процесса (может быть меньше у последнего, есл ипоровну не делится)
s1=1+myid*each                             !процесс считает от этой позиции...
s2=s1+each-1;s2=min(n,s2)                  !...до этой позиции
do s=s1,s2                                 !цикл по вершинам, перебирает правые части системы линейных уравнений. Эти правые части - орты из нулей с 1 в одной позиции
    b(:) = 0.; b(s) = 1.                   !нули везде кроме одной позиции s
    x=0.                                   !начальное предположение
    bb=b                                   !сохраним правую часть
    call pgmres(n, Krylov, b, x, vv, epsss, maxiter, 1, a, ja, ia, alu, jlu, ju, error) !вызов решателя (матрицу уже подготовили, каждый проц для себя)
    b=bb                                   !восстанавливаем правую часть
    select case(error)                     !анализируем код ошибки
    case(-1,0)                             !если все ок (или начальное приближение и есть решение)
        if(verbose) write(65,*) s,":",x(:) !запишем решение системы в файл отдельной строкой, если включен вывод
        if(check) then                     !проверяем решение системы?
            bb = 0.                        !проверка ответа. Сюда положим A*x
            do i=1,n                       !цикл по строкам
                do j = ia(i),ia(i+1)-1     !цикл по строке
                    bb(i) = bb(i) + a(j) * x(ja(j)) !собираем произведения ненулевых элементов строки и соответствующих элементов решения (№ элемента в столбце x равен номеру столбца элемента в строке)
                enddo
            enddo
            resid = maxval(abs(b-bb))      !максимум модуля разности найденного bb=A*x и правой части b (в идеале должны совпадать)
            if(resid>maxres) then          !если полученное расхожение претендует на максимальность,
                maxres  = resid            !запоминаем максимальное расхождение
                smaxres = s                ! и где оно встретилось (в параллельной версии игнорируется)
            endif
        endif
        !расчет центральности
        bb = 0.                            !здесь будут x^s(v), они же токи по вершинам i при условии, что источник в вершине s
        do i=1,n                           !цикл по строкам
            do j = ia(i),ia(i+1)-1         !цикл по строке
                bb(i) = bb(i) + abs(x(i) - x(ja(j)))*w(j) !накапливаем |phi_i^s-phi_j^s|w_{i,j}: (6) в автореферате
            enddo
            if(i.eq.s) bb(i)=bb(i)+1.      !прибавим b_s(v): см (6) в автореферате
            bb(i)=bb(i)/2.                 !делим пополам: см (6)  в автореферате
            CF(i) = CF(i)+bb(i)/n          !см (7)  в автореферате; копим для каждой вершины ток при источнике в разных вершинах s
        enddo
        print'("CPU",1X,I3,1X,":",1X,I7,1X,"out of",1X,I7,1X,"done",1X,F6.2,"% completed")', myid,s-s1+1,s2-s1+1,real(s-s1+1)/real(s2-s1+1)*100. !доклад
    case(1)                                !если сходимость не удалась
        print*, "Convergence not achieved in itmax iterations." !сообщить
        call bye                           !и умереть
    case default                           !если ХЗ что за ошибка((
        print*, "Unknown error:", error    !сообщить
        call bye                           !и умереть
    end select
enddo
if(verbose) close(65)                      !файл больше не нужен
if(myid.eq.0) allocate(CF0(n))             !создать массив чтобы собрать туда куски центральностей (их надо сложить)
call mpi_barrier(comm,mpierr)              !подождем всех (на всякий случай)
call mpi_reduce(CF,CF0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,mpierr) !собрать и сложить результаты
if(check) call mpi_reduce(maxres,maxres0,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm,mpierr) !найти максимальное отклонение
if(myid .eq. 0) then                       !отчитывается один процесс
        open(666,file=trim(fres)//'.raw')        !откроем файл для ответа в сыром виде (а то мало ли что пойдет не так)
        write(666,*) CF0(:)                !вывести центральности.
        close(666)                         !файл больше не нужен
        open(66,file=fres)                 !откроем файл для ответа
        open(99,file=fgraph)               !откроем файл со списком вершин графа: 2 колонки - номер вершины и ее метка. Делается из 3-колоночного описания графа (вершина-вершина-вес ребра) однострочником ** ниже.
        write(66,*) "===Electrical centralities===" !просто заголовок
        do i=1,n                           !цикл по вершинам графа
                read(99,*) str             !извлечь из списка вершин графа строку - номер вершины и ее метка
                write(66,'(A9,1X,F25.17)') str(2), CF0(i) !вывести метку и центральность. Сравните результаты с тестовым ответом однострочником *** 
        enddo
        close(66);close(99)                !закрыть файлы
        if(check) print*, "Maximal discrepancy",maxres0 !максимальная ошибка в решении СЛАУ - для контроля
        print*, "DONE!"                    !все хорошо
endif
call bye                                   !все убираем и уходим

CONTAINS                                   !внутренние процедуры в глобальном пространстве имен

        subroutine bye                     !уйти чисто
                deallocate(a, ja, ia, w, CF) !освободим память
                deallocate(b, bb, x)
                deallocate(vv, alu, ju, jlu) 
                deallocate(wa, jw, D)
                if(myid.eq.0) deallocate(CF0)          
                call mpi_finalize(mpierr)  !завершим MPI
                stop                       !прервем программу
        end subroutine bye

END

        !*)  3col2pack.pl
        !**) perl -n -a -E 'if(not defined $N{$F[0]}){$N{$F[0]}=$n++; say "$n $F[0]"}' EdgesOfGraph.txt > Edges.dat
        !***) perl -n -a -E 'next if /===/; push @{$h{$F[0]}}, $F[1] }{ for(sort keys %h) {$d=abs($h{$_}[0]-$h{$_}[1]);if($d>$max){$max=$d;$m=$_};say"$_ @{$h{$_}} $d"}say STDERR "max diff is $max at node $m"' res.dat centralities.txt > compare.txt

