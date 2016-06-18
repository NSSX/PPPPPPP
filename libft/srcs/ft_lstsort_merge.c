/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_lstsort_merge.c                                 :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/21 18:58:27 by igomez            #+#    #+#             */
/*   Updated: 2014/11/24 18:29:04 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

static t_list	*lst_merge(t_list *lst1, t_list *lst2,
		int (*cmp)(void *, void *))
{
	t_list *res;
	t_list *tmp;

	res = (cmp(lst1->content, lst2->content) > 0) ? lst2 : lst1;
	if ((tmp = res) == lst1)
		lst1 = lst1->next;
	else
		lst2 = lst2->next;
	while (lst1 && lst2)
	{
		if (cmp(lst1->content, lst2->content) > 0)
		{
			tmp->next = lst2;
			lst2 = lst2->next;
		}
		else
		{
			tmp->next = lst1;
			lst1 = lst1->next;
		}
		tmp = tmp->next;
	}
	tmp->next = (lst1) ? lst1 : lst2;
	return (res);
}

static void		lst_sort(t_list **ptr_lst, int (*cmp)(void *, void *))
{
	t_list	*lst2;
	t_list	*slow;
	t_list	*fast;

	lst2 = *ptr_lst;
	if (!lst2 || !(lst2->next))
		return ;
	fast = lst2;
	while (fast && fast->next)
	{
		slow = lst2;
		fast = fast->next->next;
		lst2 = lst2->next;
	}
	slow->next = NULL;
	lst_sort(ptr_lst, cmp);
	lst_sort(&lst2, cmp);
	*ptr_lst = lst_merge(*ptr_lst, lst2, cmp);
}

void			ft_lstsort_merge(t_list **lst_ptr,
					int (*cmp)(void *, void *))
{
	if (!lst_ptr || !(*lst_ptr) || !cmp)
		return ;
	lst_sort(lst_ptr, cmp);
}
